module VoronoiRemeshModule

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use TracerAndVorticityDistributionModule
use ParticlesEdgesPanelsModule
use RefineRemeshModule
use VoronoiPanelsModule

implicit none

private
public RemeshVoronoi
!public SINGLE_GAUSSIAN_VORTEX, RH4_WAVE, MULTIPLE_VORTICES, RH2_WAVE, &
!	   STRATOSPHERE_MODEL, JET, ZONAL_MEAN

!integer(kint), parameter :: SINGLE_GAUSSIAN_VORTEX = 41, &
!							RH4_WAVE = 42,&
!							MULTIPLE_VORTICES = 43,&
!							RH2_WAVE = 44, &
!							STRATOSPHERE_MODEL = 45, &
!							JET = 46, &
!							ZONAL_MEAN = 47

integer(kint), parameter :: logLevel=DEBUG_LOGGING_LEVEL
type(Logger) :: log
character(len=28) :: logKey='VoronoiRemesh : '
character(len=256) :: logString
logical(klog), save :: logInit = .False.

contains

subroutine RemeshVoronoi(aPanels,initNest,AMR,amrTol1,amrTol2,procRank,problemID,&
	theta0,beta,perturbAmp,perturbWaveNum, time)
	type(VorPanels), pointer, intent(inout) :: aPanels
	integer(kint), intent(in) :: initNest, AMR, procRank, problemID
	real(kreal), intent(in) :: amrTol1, amrTol2
	real(kreal), intent(in), optional :: theta0, beta, perturbAmp, time
	integer(kint), intent(in), optional :: perturbWaveNum
	! SSRFPACK Variables
	real(kreal), allocatable :: gradx0(:,:), grady0(:,:), gradz0(:,:)
	real(kreal), allocatable :: sigmax0(:), sigmay0(:), sigmaz0(:)
	real(kreal) :: sigmaTol, dSigX, dSigY, dSigZ
	integer(kint) :: errCode, gradFlag, sigmaFlag
	! Remesh variables
	type(VorPanels), pointer :: newPanels, tempPanels
	type(Particles), pointer :: triParticles
	type(Edges), pointer :: triEdges
	type(Panels), pointer :: triPanels
	logical(klog), allocatable :: refineFlag(:)
	integer(kint) :: nTracer, startPanel
	real(kreal) :: newLat, newLon, norm
	integer(kint) :: i, j, k, n
	integer(kint) :: amrCount1, amrCount2, amrLoopCounter, refineCount, startIndex
	integer(kint) :: nOldPanels, nOldPanels2, nOldParticles, spaceLeft, vertIndices(3)
	real(kreal) :: maxx0(3), minx0(3), lagVar
	logical(klog) :: maxRefinement
	
	maxRefinement = .False.
	if ( .NOT. logInit) call InitLogger(log)
	
	call LogMessage(log,DEBUG_LOGGING_LEVEL,logKey,'Entering RemeshVoronoi.')
	
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 1 : Setup existing mesh as source for data interpolation  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	
	n = aPanels%N
	! Allocate variables for SSRFPACK
	allocate(gradx0(3,n))
	allocate(grady0(3,n))
	allocate(gradz0(3,n))
	allocate(sigmax0(6*n-12))
	allocate(sigmay0(6*n-12))
	allocate(sigmaz0(6*n-12))
	
	! Estimate the gradients at each generator using Delaunay triangulation and SSRFPACK
	do j=1,n
		call GRADL(n,j,aPanels%x,aPanels%y,aPanels%z,aPanels%x0, &
				   aPanels%list,aPanels%lptr,aPanels%lend,gradx0(:,j),errCode)
		call GRADL(n,j,aPanels%x,aPanels%y,aPanels%z,aPanels%y0, &
				   aPanels%list,aPanels%lptr,aPanels%lend,grady0(:,j),errCode)
		call GRADL(n,j,aPanels%x,aPanels%y,aPanels%z,aPanels%z0, &
				   aPanels%list,aPanels%lptr,aPanels%lend,gradz0(:,j),errCode)		   		   
	enddo
	gradFlag = 1
	
	! Compute the shape-preserving tension factors
	sigmaTol = 0.01_kreal
	sigmax0 = 0.0_kreal
	sigmay0 = 0.0_kreal
	sigmaz0 = 0.0_kreal
	call GETSIG(n,aPanels%x,aPanels%y,aPanels%z,aPanels%x0, & 
				aPanels%list,aPanels%lptr,aPanels%lend, gradx0, sigmaTol, sigmax0,dSigX,errCode)
	call GETSIG(n,aPanels%x,aPanels%y,aPanels%z,aPanels%y0, & 
				aPanels%list,aPanels%lptr,aPanels%lend, grady0, sigmaTol, sigmay0,dSigY,errCode)	
	call GETSIG(n,aPanels%x,aPanels%y,aPanels%z,aPanels%z0, & 
				aPanels%list,aPanels%lptr,aPanels%lend, gradz0, sigmaTol, sigmaz0,dSigZ,errCode)
	sigmaFlag = 1				
	
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 2 : Build new mesh and use as destination for interpolation  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	
	nTracer = GetNTRacer(aPanels)
	if ( AMR == 0 ) then
		call New(newPanels,initNest,AMR,nTracer)
		startPanel = 1
		! Interpolate x0, y0, z0 from old mesh to new mesh
		do j=1, n
			newLat = Latitude(newPanels%x(j),newPanels%y(j),newPanels%z(j))
			newLon = Longitude(newPanels%x(j), newPanels%y(j), newPanels%z(j))
			call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z, aPanels%x0, &
						aPanels%list,aPanels%lptr,aPanels%lend, sigmaFlag, sigmax0, &
						gradFlag, gradx0, startPanel, newPanels%x0(j), errCode)
			call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z, aPanels%y0, &
						aPanels%list,aPanels%lptr,aPanels%lend, sigmaFlag, sigmay0, &
						gradFlag, grady0, startPanel, newPanels%y0(j), errCode)	
			call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z, aPanels%z0, &
						aPanels%list,aPanels%lptr,aPanels%lend, sigmaFlag, sigmaz0, &
						gradFlag, gradz0, startPanel, newPanels%z0(j), errCode)					
		enddo
		do j=1,n
			norm = sqrt( newPanels%x0(j)**2 + newPanels%y0(j)**2 + newPanels%z0(j)**2)
			newPanels%x0(j) = newPanels%x0(j)/norm
			newPanels%y0(j) = newPanels%y0(j)/norm
			newPanels%z0(j) = newPanels%z0(j)/norm
		enddo
		!!!!!
		! SET VORTICITY
		!!!!!	
		! Set absolute vorticity on new grid
			if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
				do j=1,n
					newPanels%absVort(j) = GaussianVortexX(newPanels%x0(j),newPanels%y0(j),newPanels%z0(j)) - &
						gaussConst + 2.0_kreal*Omega*newPanels%z0(j)
				enddo
			elseif (problemID == RH4_WAVE) then
				do j=1,n
					newPanels%absVort(j) = HaurwitzStationary4RelVort(newPanels%x0(j),newPanels%y0(j),newPanels%z0(j)) + &
						2.0_kreal*Omega*newPanels%z0(j)
				enddo
			elseif (problemID == MULTIPLE_VORTICES) then
				do j=1,n
					newPanels%absVort(j) = ManyVortsX(newPanels%x0(j),newPanels%y0(j),newPanels%z0(j)) - &
						gaussConst + 2.0_kreal*Omega*newPanels%z0(j)
				enddo
			elseif ( problemID == RH2_WAVE ) then
				do j=1,n
					newPanels%absVort(j) = RH2WaveX(newPanels%x0(j),newPanels%y0(j),newPanels%z0(j)) + &
										   2.0_kreal*Omega*newPanels%z0(j)
				enddo
			elseif ( problemID == STRATOSPHERE_MODEL) then
			
			elseif ( problemID == JET ) then	
				do j=1,n
					newPanels%absVort(j) = JetRelVortX(newPanels%x0(j),newPanels%y0(j),newPanels%z0(j),&
						theta0,beta,perturbAmp,perturbWaveNum) + &
						2.0_kreal*Omega*newPanels%z0(j) - gaussConst
				enddo
			elseif ( problemID == ZONAL_MEAN) then
				do j=1,n
					newPanels%absVort(j) = ZonalMeanX(newPanels%x0(j),newPanels%y0(j),newPanels%z0(j)) + &
										   2.0_kreal*Omega*newPanels%z0(j)
				enddo
			endif
		! Set relative vorticity on new grid	
			do j=1,n
				newPanels%relVort(j) = newPanels%absVort(j) - 2.0_kreal*Omega*newPanels%z(j)
			enddo	
	else! VORONOI AMR
		 call New(triParticles,triEdges,triPanels,TRI_PANEL,initNest,AMR,nTracer,BVE_SOLVER)
		 allocate(refineFlag(triPanels%N_Max))
		 startPanel=1
		 refineFlag = .FALSE.
		 do j=1,triPanels%N
		 	if ( .NOT. triPanels%hasChildren(j) ) then
		 		newLat = Latitude(triPanels%x(:,j))
		 		newLon = Longitude(triPanels%x(:,j))
		 		call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%x0,&
		 			aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmax0,&
		 			gradFlag,gradx0,startPanel,triPanels%x0(1,j),errCode)
		 		call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%y0,&
		 			aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmay0,&
		 			gradFlag,grady0,startPanel,triPanels%x0(2,j),errCode)	
		 		call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%z0,&
		 			aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmaz0,&
		 			gradFlag,gradz0,startPanel,triPanels%x0(3,j),errCode)
		 		triPanels%x0(:,j) = triPanels%x0(:,j)/sqrt(sum(triPanels%x0(:,j)*triPanels%x0(:,j)))	
		 	endif
		 enddo
		 startPanel=1
		 do j=1,triParticles%N
		 	newLat = Latitude(triParticles%x(:,j))
		 	newLon = Longitude(triParticles%x(:,j))
		 	call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%x0,&
		 		aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmax0,&
		 		gradFlag,gradx0,startPanel,triParticles%x0(1,j),errCode)
		 	call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%y0,&
		 		aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmay0,&
		 		gradFlag,grady0,startPanel,triParticles%x0(2,j),errCode)
		 	call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%z0,&
		 		aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmaz0,&
		 		gradFlag,gradz0,startPanel,triParticles%x0(3,j),errCode)
		 	triParticles%x0(:,j) = triParticles%x0(:,j)/sqrt(sum(triParticles%x0(:,j)*triParticles%x0(:,j)))	
		 enddo
		 if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
			do k=1,triParticles%N
				triParticles%absVort(k) = GaussianVortexX(triParticles%x0(:,k)) - gaussConst + &
					2.0_kreal*Omega*triParticles%x0(3,k)
			enddo
			do k=1,triPanels%N
				if ( .NOT. triPanels%hasChildren(k) ) then
					triPanels%absVort(k) = GaussianVortexX(triPanels%x0(:,k)) - gaussConst + &
						2.0_kreal*Omega*triPanels%x0(3,k)
				endif
			enddo
		elseif (problemID == RH4_WAVE ) then
			do k=1,triParticles%N
				triParticles%absVort(k) = HaurwitzStationary4RelVort(triParticles%x0(:,k)) + &
					 2.0_kreal*Omega*triParticles%x0(3,k)
			enddo
			do k=1,triPanels%N
				if ( .NOT. triPanels%hasCHildren(k) ) then
					triPanels%absVort(k) = HaurwitzStationary4RelVort(triPanels%x0(:,k)) + &
						2.0_kreal*Omega*triPanels%x0(3,k)
				endif
			enddo					
		elseif ( problemID == MULTIPLE_VORTICES ) then
			do k=1,triParticles%N
				triParticles%absVort(k) = ManyVortsX(triParticles%x0(:,k)) - gaussConst + &
					2.0_kreal*Omega*triParticles%x0(3,k)
			enddo
			do k=1,triPanels%N
				if ( .NOT. triPanels%hasChildren(k) ) then
					triPanels%absVort(k) = ManyVortsX(triPanels%x0(:,k)) - gaussConst + &
						2.0_kreal*Omega*triPanels%x0(3,k)
				endif		
			enddo
		elseif ( problemID == RH2_WAVE) then
			do k=1,triParticles%N
				triParticles%absVort(k) = RH2WaveX(triParticles%x0(:,k)) + &
						2.0_kreal*Omega*triParticles%x0(3,k)
			enddo
			do k=1,triPanels%N
				if ( .NOT. triPanels%hasChildren(k) ) then
					triPanels%absVort(k) = RH2WaveX(triPanels%x0(:,k)) + &
						2.0_kreal*Omega*triPanels%x0(3,k)
				endif
			enddo
		elseif ( problemID == STRATOSPHERE_MODEL ) then
			if ( .NOT. present(time) ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Remesh Stratosphere ERROR : time argument missing.')
				return
			endif
			do k=1,triParticles%N
				triParticles%absVort(k) = StratosphereRelVortX(triParticles%x0(:,k)) - gaussConst &
					+ 2.0_kreal*Omega*triParticles%x0(3,k)
			enddo
			do k=1,triPanels%N
				if ( .NOT. triPanels%hasChildren(k)) then
					triPanels%absVort(k) = StratosphereRelVortX(triPanels%x0(:,k)) - gaussConst &
						+ 2.0_kreal*Omega*triPanels%x0(3,k)
				 endif
			enddo
		elseif ( problemID == JET ) then
			if ( (( .NOT. present(perturbWaveNum) ) .OR. (.NOT. present(beta))) .OR. &
				 ( (.NOT. present(theta0)) .OR. (.NOT. present(perturbAmp)) ) ) then
				 call LogMessage(log,ERROR_LOGGING_LEVEL,logkey,'Remesh Jet ERROR : JetRelVortX arguments needed.')
				 return
			endif
			do k=1,triParticles%N
				triParticles%absVort(k) = JetRelVortX(triParticles%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
					2.0_kreal*Omega*triParticles%x0(3,k) - gaussConst
			enddo
			do k=1,triPanels%N
				if (.NOT. triPanels%hasChildren(k) ) then
					triPanels%absVort(k) = JetRelVortX(triPanels%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
						2.0_kreal*Omega*triPanels%x0(3,k) - gaussConst
				endif
			enddo
		elseif ( problemID == ZONAL_MEAN) then
			do k=1,triParticles%N
				triParticles%absVort(k) = ZonalMeanX(triParticles%x0(:,k)) + &
					2.0_kreal*Omega*triParticles%x0(3,k)
			enddo
			do k=1,triPanels%N
				if ( .NOT. triPanels%hasChildren(k) ) then
					triPanels%absVort(k) = ZonalMeanX(triPanels%x0(:,k)) + &
						2.0_kreal*Omega*triPanels%x0(3,k)
				endif
			enddo
		elseif ( problemID == GAUSSIAN_HILLS) then
			do k=1,triParticles%N
				triParticles%tracer(k,1) = GaussianHillsX(triParticles%x0(:,k))
			enddo
			do k=1,triPanels%N
				if (.NOT. triPanels%hasChildren(k) ) then
					triPanels%tracer(k,1) = GaussianHillsX(triPanels%x0(:,k))
				endif
			enddo
		elseif ( problemID == COSINE_BELLS) then
			do k=1,triParticles%N
				triParticles%tracer(k,1) = CosineBellsX(triParticles%x0(:,k))
			enddo
			do k=1,triPanels%N
				if (.NOT. triPanels%hasChildren(k) ) then
					triPanels%tracer(k,1) = CosineBellsX(triPanels%x0(:,k))
				endif
			enddo
		elseif ( problemID == SLOTTED_CYLINDERS) then		
			do k=1,triParticles%N
				triParticles%tracer(k,1) = SlottedCylindersX(triParticles%x0(:,k))
			enddo
			do k=1,triPanels%N
				if (.NOT. triPanels%hasChildren(k) ) then
					triPanels%tracer(k,1) = SlottedCylindersX(triPanels%x0(:,k))
				endif
			enddo
		endif
		if ( problemID == STRATOSPHERE_MODEL ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'StratModel not implemented yet.')
			return
		else
			triParticles%relVort(1:triParticles%N) = triParticles%absVort(1:triParticles%N) &
				- 2.0_kreal*Omega*triParticles%x(3,1:triParticles%N)
			do j=1,triPanels%N
				if (.NOT. triPanels%hasChildren(j) ) then
					triPanels%relVort(j) = triPanels%absVort(j) - 2.0_kreal*Omega*triPanels%x(3,j)
				endif
			enddo
		endif
		! UNIFORM Triangular mesh ready
		amrCount1 = 0
		amrCount2 = 0
		do j=1,triPanels%N
			if (.NOT. triPanels%hasChildren(j) ) then
				if ( (problemID>=GAUSSIAN_HILLS) .AND. (problemID<=SLOTTED_CYLINDERS)) then
					call LogMessage(log,WARNING_LOGGING_LEVEL,logkey,'Advection AMR not implemented yet.')
					exit
				else
!!!!!!
! AMR CRITERIA
!!!!!!				amr criterion 1
					if ( abs(triPanels%relVort(j))*triPanels%area(j) > amrTol1) then
						refineFlag(j) = .TRUE.
						amrCount1 = amrCount1 + 1
					else
						maxX0 = triPanels%x0(:,j)
						minX0 = triPanels%x0(:,j)
						do k=1,3
							if ( triParticles%x0(1,triPanels%vertices(k,j)) < minX0(1) ) then
								minX0(1) = triParticles%x0(1,triPanels%vertices(k,j))
							endif
							if ( triParticles%x0(1,triPanels%vertices(k,j)) > maxX0(1) ) then
								maxX0(1) = triParticles%x0(1,triPanels%vertices(k,j))
							endif
							if ( triParticles%x0(2,triPanels%vertices(k,j)) < minx0(2)) then
								minX0(2) = triParticles%x0(2,triPanels%vertices(k,j))
							endif
							if ( triParticles%x0(2,triPanels%vertices(k,j)) > maxx0(2)) then
								maxX0(2) = triParticles%x0(2,triPanels%vertices(k,j))
							endif
							if ( triParticles%x0(3,triPanels%vertices(k,j)) < minx0(3)) then
								minX0(3) = triParticles%x0(3,triPanels%vertices(k,j))
							endif
							if ( triParticles%x0(3,triPanels%vertices(k,j)) > maxx0(3)) then
								maxX0(3) = triParticles%x0(3,triPanels%vertices(k,j))
							endif
						enddo
						lagVar = sum(maxx0-minx0)
						if ( lagVar > amrTol2) then
							refineFlag(j) = .TRUE.
							amrCount2 = amrCount2+1
						endif
					endif!amr criteria		
				endif!bve problem
			endif!low-level panel
		enddo!initial triPanels loop
		refineCount = count(refineFlag)
		startIndex = 1
		amrLoopCounter = 1
		startPanel=1
		do while (refineCount > 0) 
			spaceLeft = triPanels%N_Max - triPanels%N
			if (spaceLeft/4 > refineCount ) then
				if ( procRank <= 0 ) then
					write(logString,'(A,I4,A,I8,A)') " AMR Loop ",amrLoopCounter,&
						": refining ",refineCount," panels."
					call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),logString)
				endif
				! loop over grid, refine flagged panels
				nOldPanels = triPanels%N
				do j=startIndex,nOldPanels
					if ( refineFlag(j) ) then
						! Divide panel
						nOldParticles = triParticles%N
						nOldPanels2 = triPanels%N
						call DividePanel(triParticles,triEdges,triPanels,j)
						refineFlag(j) = .FALSE.
						triPanels%absVort(j) = 0.0_kreal
						triPanels%relVort(j) = 0.0_kreal
						triPanels%area(j) = 0.0_kreal
						do k=nOldParticles+1,triParticles%N
							newLat = Latitude(triParticles%x(:,k))
							newLon = Longitude(triParticles%x(:,k))
							call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%x0,&
								aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmax0,&
								gradFlag,gradx0,startPanel,triParticles%x0(1,k),errCode)
							call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%y0,&
								aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmay0,&
								gradFlag,grady0,startPanel,triParticles%x0(2,k),errCode)
							call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%z0,&
								aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmaz0,&
								gradFlag,gradz0,startPanel,triParticles%x0(3,k),errCode)
							triParticles%x0(:,k) = triParticles%x0(:,k)/sqrt(sum(triParticles%x0(:,k)*triParticles%x0(:,k)))
!!!!!
! SET VORTICITY					Set vorticity on new particles
!!!!!							
							if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
								triParticles%absVort(k) = GaussianVortexX(triParticles%x0(:,k)) - &
									gaussConst + 2.0_kreal*Omega*triParticles%x0(3,k)
							elseif (problemID == RH4_WAVE ) then
								triParticles%absVort(k) =  HaurwitzStationary4RelVort(triParticles%x0(:,k)) +&
											 2.0_kreal*Omega*triParticles%x0(3,k)																		
							elseif ( problemID == MULTIPLE_VORTICES ) then
								triParticles%absVort(k) = ManyVortsX(triParticles%x0(:,k)) - &
									gaussConst + 2.0_kreal*Omega*triParticles%x0(3,k)
							elseif ( problemID == RH2_Wave ) then
								triParticles%absVort(k) = RH2WaveX(triParticles%x0(:,k)) + &
									2.0_kreal*Omega*triParticles%x0(3,k)
							elseif ( problemID == STRATOSPHERE_MODEL) then
								triParticles%absVort(k) = StratosphereRelVortX(triParticles%x0(:,k)) - gaussConst &
									+ 2.0_kreal*Omega*triParticles%x0(3,k)
							elseif ( problemID == JET) then
								triParticles%absVort(k) = JetRelVortX(triParticles%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
									2.0_kreal*Omega*triParticles%x0(3,k) - gaussConst
							elseif ( problemID == ZONAL_MEAN ) then
								triParticles%absVort(k) = ZonalMeanX(triParticles%x0(:,k)) + &
									2.0_kreal*Omega*triParticles%x0(3,k)
							elseif ( problemID == GAUSSIAN_HILLS) then
								triParticles%tracer(k,1) = GaussianHillsX(triParticles%x0(:,k))
							elseif ( problemID == COSINE_BELLS) then
								triParticles%tracer(k,1) = CosineBellsX(triParticles%x0(:,k))
							elseif ( problemID == SLOTTED_CYLINDERS) then										
								triParticles%tracer(k,1) = SlottedCylindersX(triParticles%x0(:,k))
							endif
							
							if ( problemID == STRATOSPHERE_MODEL) then
								triParticles%relVort(k) = triParticles%absVort(k) - &
									2.0_kreal*Omega*triParticles%x(3,k) &
									- Juckes_Forcing(triParticles%x(:,k),time)
							else
								triParticles%relVort(k) = triParticles%absVort(k) - &
										2.0_kreal*Omega*triParticles%x(3,k)
							endif
						enddo! new particles
						! find x0 of new panels
						do k=nOldPanels2+1,triPanels%N
							newLat = Latitude(triPanels%x(:,k))
							newLon = Longitude(triPanels%x(:,k))
							call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%x0,&
								aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmax0,&
								gradFlag,gradx0,startPanel,triPanels%x0(1,k),errCode)
							call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%y0,&
								aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmay0,&
								gradFlag,grady0,startPanel,triPanels%x0(2,k),errCode)	
							call INTRC1(n,newLat,newLon,aPanels%x,aPanels%y,aPanels%z,aPanels%z0,&
								aPanels%list,aPanels%lptr,aPanels%lend,sigmaFlag,sigmaz0,&
								gradFlag,gradz0,startPanel,triPanels%x0(3,k),errCode)
							triPanels%x0(:,k) = triPanels%x0(:,k)/sqrt(sum(triPanels%x0(:,k)*triPanels%x0(:,k)))
	!!!!!
	! SET VORTICITY				Set vorticity on new panels
	!!!!!						
							if ( problemID == SINGLE_GAUSSIAN_VORTEX) then
								triPanels%absVort(k) = GaussianVortexX(triPanels%x0(:,k)) - &
									gaussConst + 2.0_kreal*Omega*triPanels%x0(3,k)
							elseif (problemID == RH4_WAVE ) then		
								triPanels%absVort(k) = HaurwitzStationary4RelVort(triPanels%x0(:,k)) + &
									2.0_kreal*Omega*triPanels%x0(3,k)
							elseif ( problemID == MULTIPLE_VORTICES ) then
								triPanels%absVort(k) = ManyVortsX(triPanels%x0(:,k)) - gaussConst + &
									2.0_kreal*Omega*triPanels%x0(3,k)
							elseif ( problemID == RH2_WAVE ) then
								triPanels%absVort(k) = RH2WaveX(triPanels%x0(:,k)) + &
									2.0_kreal*Omega*triPanels%x0(3,k)
							elseif ( problemID == STRATOSPHERE_MODEL ) then
								triPanels%absVort(k) = StratosphereRelVortX(triPanels%x0(:,k)) - gaussConst &
									+ 2.0_kreal*Omega*triPanels%x0(3,k)
							elseif ( problemID == JET ) then
								triPanels%absVort(k) = JetRelVortX(triPanels%x0(:,k),theta0,beta,perturbAmp,perturbWaveNum) + &
									2.0_kreal*Omega*triPanels%x0(3,k) - gaussConst
							elseif ( problemID == ZONAL_MEAN ) then
								triPanels%absVort(k) = ZonalMeanX(triPanels%x0(:,k)) + &
									2.0_kreal*Omega*triPanels%x0(3,k)
							elseif ( problemID == GAUSSIAN_HILLS) then
								triPanels%tracer(k,1) = GaussianHillsX(triPanels%x0(:,k))
							elseif ( problemID == COSINE_BELLS) then
								triPanels%tracer(k,1) = CosineBellsX(triPanels%x0(:,k))
							elseif ( problemID == SLOTTED_CYLINDERS) then
								triPanels%tracer(k,1) = SlottedCylindersX(triPanels%x0(:,k))
							endif
							
							if ( problemID == STRATOSPHERE_MODEL ) then
								triPanels%relVort(k) = triPanels%absVort(k) - &
										2.0_kreal*Omega*triPanels%x(3,k) &
										- Juckes_Forcing(triPanels%x(:,k),time)							
							else		
								triPanels%relVort(k) = triPanels%absVort(k) - &
										2.0_kreal*Omega*triPanels%x(3,k)
							endif
	!!!!!
	! AMR CRITERIA				Check refinement criteria on new panels
	!!!!!					    
							if ( triPanels%nest(j) <= initNest + MAX_REFINEMENT) then
								if ( (problemID >= GAUSSIAN_HILLS) .AND. ( problemID <= SLOTTED_CYLINDERS) ) then
	! ADVECTION			
									if ( abs(triPanels%tracer(k,1))*triPanels%area(k) > amrTol1 ) then
										refineFlag(k) = .TRUE.
										amrCount1 = amrCount1 + 1
									endif
								else
	! BVE								
									if ( abs(triPanels%relVort(k))*triPanels%area(k) > amrTol1 ) then
										refineFlag(k) = .TRUE.
										amrCount1 = amrCount1 + 1
									else
										minx0 = triPanels%x0(:,k)
										maxx0 = triPanels%x0(:,k)
										vertIndices = triPanels%vertices(:,k)
										do i=1,3
											if ( triParticles%x0(1,vertIndices(i)) < minx0(1) ) then
												minx0(1) = triParticles%x0(1,vertIndices(i))
											endif
											if ( triParticles%x0(2,vertIndices(i)) < minx0(2) ) then
												minx0(2) = triParticles%x0(2,vertIndices(i))
											endif
											if ( triParticles%x0(3,vertIndices(i)) < minx0(3) ) then
												minx0(3) = triParticles%x0(3,vertIndices(i))
											endif
											if ( triParticles%x0(1,vertIndices(i)) > maxx0(1) ) then
												maxx0(1) = triParticles%x0(1,vertIndices(i))
											endif
											if ( triParticles%x0(2,vertIndices(i)) > maxx0(2) ) then
												maxx0(2) = triParticles%x0(2,vertIndices(i))
											endif
											if ( triParticles%x0(3,vertIndices(i)) > maxx0(3) ) then
												maxx0(3) = triParticles%x0(3,vertIndices(i))
											endif
										enddo
								
										lagVar = sum(maxx0-minx0)
										if ( lagVar > amrTol2 ) then
											refineFlag(k) = .TRUE.
											amrCount2 = amrCount2 + 1
										endif
									endif!amr 2
								endif
							else
								if ( .NOT. maxRefinement) then
									if ( procRank == 0 ) call LogMessage(log,WARNING_LOGGING_LEVEL,logKey,'Refinement WARNING : maximum refinement reached.')
									maxRefinement = .TRUE.
								endif
							endif
						enddo!newPanels
					endif!flagged panel
				enddo!whole panels loop
				refineCount = count(refineFlag)
				amrLoopCounter = amrLoopCounter + 1
				startIndex = nOldPanels
			else
				if ( procRank <= 0 ) call LogMessage(log,WARNING_LOGGING_LEVEL,trim(logKey),' WARNING: not enough memory for AMR.')
				exit
			endif ! spaceleft
		enddo!while refineCount > 0
		! Log statistics
		if ( procRank <= 0) then
			write(logString,'(A,I8,A)') " AMR circulation criterion triggered ",amrCount1," times."
			call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),trim(logString))
			write(logString,'(A,I8,A)') " AMR Lag. coord. variation criterion triggered ",amrCount2," times."
			call LogMessage(log,TRACE_LOGGING_LEVEL,trim(logKey),trim(logString))
		endif
		
		call New(newPanels,triParticles%N,nTracer,triParticles)
		call Delete(triParticles,triEdges,triPanels)
		deallocate(refineFlag)
endif
			
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Part 3 : Replace old grid with new grid  !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	tempPanels => aPanels
	
	aPanels => newPanels
	
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
	!	Delete old grid & Clean up	    !
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

	call Delete(tempPanels)
	nullify(newPanels)
						
	deallocate(gradx0)
	deallocate(grady0)
	deallocate(gradz0)
	deallocate(sigmax0)
	deallocate(sigmay0)
	deallocate(sigmaz0)
end subroutine


subroutine NCLOutputLLPlot(self,fileroot, degreeSpacing)
	type(VorPanels), intent(in) :: self
	character(len=*), intent(in) :: fileroot
	real(kreal), intent(in) :: degreeSpacing
	character(len=256) :: filename
	character(len=28) :: datastring
	integer(kint) :: nTracer
	real(kreal), allocatable :: lats(:), lons(:), interp(:,:), sigma(:), grad(:,:)
	real(kreal) :: dLambda, sigmaTol, dSig
	integer(kint) :: nLat, nLon, i, j, k, sigmaFlag, gradientFlag, errCode
	integer(kint), parameter :: writeUnit = 11
	integer(kint) :: writeStat
	
	if ( associated(self%tracer) ) then
		nTracer = size(self%tracer,2)
	else
		nTracer = 0
	endif
	
	nLat = floor(180.0_kreal/degreeSpacing) + 1
	nLon = floor(360.0_kreal/degreeSpacing)
	dLambda = degreeSpacing*PI/180.0_kreal
	
	allocate(lats(nLat))
	allocate(lons(nLon))
	allocate(interp(nLat,nLon))
	allocate(sigma(6*self%N-12))
	allocate(grad(3,self%N))
	
	
	do j=1,nLat
		lats(j) = -PI/2.0_kreal + (j-1)*dLambda
	enddo
	do j=1,nLon
		lons(j) = (j-1)*dLambda
	enddo
	interp = 0.0_kreal
	
! Interpolate Absolute vorticity to LL grid
	write(filename,'(A,A)') trim(fileRoot),'_absVort.dat'
	
	! Estimate gradients
	do j=1,self%N
		call GRADL(self%N,j,self%x,self%y,self%z,self%absVort,&
				   self%list,self%lptr,self%lend,grad(:,j),errCode)
	enddo
	gradientFlag = 1

	! Find tension factors
	sigma = 0.0_kreal
	sigmaTol = 0.1_kreal
	call GETSIG(self%N,self%x,self%y,self%z,self%absVort, &
			    self%list,self%lptr,self%lend,grad,sigmaTol,sigma,dSig,errCode)
	sigmaFlag = 1
	
	call UNIF(self%N,self%x,self%y,self%z,self%absVort,&
			  self%list,self%lptr,self%lend,sigmaFlag,sigma,&
			  nLat,nLat,nLon,lats,lons,gradientFlag,grad,interp,errCode)
	if ( errCode == -1 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK invalid arguments.')
		return
	elseif ( errCode == -2 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK found collinear nodes.')
		return
	elseif ( errCode == -3 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK extrapolation failed.')
		return
	endif
	
	open(unit=writeUnit,file=filename,action='WRITE',status='REPLACE',iostat=errCode)
		if (errCode /= 0 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : unable to open absVort file.')
			return
		endif
		do i=1,nLat
			do j=1,nLon-1
				write(writeUnit,'(F24.15)',advance='NO') interp(i,j)
			enddo
			write(writeUnit,'(F24.15)',advance='YES') interp(i,nLon)
		enddo
	close(writeUnit)
			  
! Relative vorticity interpolation
	write(filename,'(A,A)') trim(fileRoot),'_relVort.dat'
	do j=1,self%N
		call GRADL(self%N,j,self%x,self%y,self%z,self%relVort,&
				   self%list,self%lptr,self%lend,grad(:,j),errCode)
	enddo
	! Find tension factors
	sigma = 0.0_kreal
	sigmaTol = 0.1_kreal
	call GETSIG(self%N,self%x,self%y,self%z,self%relVort, &
			    self%list,self%lptr,self%lend,grad,sigmaTol,sigma,dSig,errCode)
	call UNIF(self%N,self%x,self%y,self%z,self%relVort,&
			 self%list,self%lptr,self%lend,sigmaFlag,sigma,&
			  nLat,nLat,nLon,lats,lons,gradientFlag,grad,interp,errCode)
	if ( errCode == -1 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK invalid arguments.')
		return
	elseif ( errCode == -2 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK found collinear nodes.')
		return
	elseif ( errCode == -3 ) then
		call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK extrapolation failed.')
		return
	endif
	open(unit=writeUnit,file=filename,action='WRITE',status='REPLACE',iostat=errCode)
		if (errCode /= 0 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : unable to open relVort file.')
			return
		endif
		do i=1,nLat
			do j=1,nLon-1
				write(writeUnit,'(F24.15)',advance='NO') interp(i,j)
			enddo
			write(writeUnit,'(F24.15)',advance='YES') interp(i,nLon)
		enddo
	close(writeUnit)

! Tracer interpolation
	do k=1,nTracer
		write(dataString,'(A,I1,A)') '_tracer',k,'.dat'
		write(filename,'(A,A)') trim(fileroot),trim(datastring)
		do j=1,self%N
			call GRADL(self%N,j,self%x,self%y,self%z,self%tracer(:,k),&
				   self%list,self%lptr,self%lend,grad(:,j),errCode)
		enddo
		sigma = 0.0_kreal
		sigmaTol = 0.1_kreal
		call GETSIG(self%N,self%x,self%y,self%z,self%tracer(:,k), &
					self%list,self%lptr,self%lend,grad,sigmaTol,sigma,dSig,errCode)
		call UNIF(self%N,self%x,self%y,self%z,self%tracer(:,k),&
			 self%list,self%lptr,self%lend,sigmaFlag,sigma,&
			  nLat,nLat,nLon,lats,lons,gradientFlag,grad,interp,errCode)
		if ( errCode == -1 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK invalid arguments.')
			return
		elseif ( errCode == -2 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK found collinear nodes.')
			return
		elseif ( errCode == -3 ) then
			call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : SSRFPACK extrapolation failed.')
			return
		endif			
		open(unit=writeUnit,file=filename,action='WRITE',status='REPLACE',iostat=errCode)
			if (errCode /= 0 ) then
				call LogMessage(log,ERROR_LOGGING_LEVEL,logKey,'NCLOutput ERROR : unable to open tracer file.')
				return
			endif
			do i=1,nLat
				do j=1,nLon-1
					write(writeUnit,'(F24.15)',advance='NO') interp(i,j)
				enddo
				write(writeUnit,'(F24.15)',advance='YES') interp(i,nLon)
			enddo
		close(writeUnit)			
	enddo

			
	deallocate(sigma)
	deallocate(grad)
	deallocate(interp)
	deallocate(lats)
	deallocate(lons)
end subroutine


subroutine InitLogger(aLog)
	type(Logger), intent(inout) :: aLog
	integer(kint) :: logUnit
	logUnit = 6
	call New(aLog,logLevel,logUnit)
	logInit = .True.
end subroutine



end module