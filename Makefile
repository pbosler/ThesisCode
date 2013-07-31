## Makefile for SphereGrids3 Project

##########################
#  	Default	Compilers  	 #
##########################
FF=ifort
FF_FLAGS=-O2 
CPP=icpc
CPP_FLAGS=-O2

##########################
#  LIBRARIES (Ferrari) 	 #
##########################
MKL_INCLUDE=-I/opt/intel/composer_xe_2011_sp1.11.344/mkl/include
MKL_LIB_DIR=-L/opt/intel/composer_xe_2011_sp1.11.344/mkl/lib
MKL_LIBS=-lmkl_rt
VTK_LIB_DIR=-L/usr/local/lib/vtk-5.10
VTK_INCLUDE=-I/usr/local/include/vtk-5.10
VTK_LIBS=-lvtkCommon -lvtkGraphics -lvtkRendering -lvtkViews -lvtkWidgets -lvtkImaging \
	-lvtkHybrid -lvtkIO -lvtkFiltering
NCL_LIB_DIR=-L/Users/pbosler/SourceCodes/ncl_ncarg-6.0/lib
NCL_LIBS=-lngmath
GFORTRAN_LIB_DIR=-L/usr/local/gfortran/lib
GFORTRAN_LIBS=-lgfortran

########################################
# 	Pattern rules and Make variables   #
########################################
BASE_OBJS=NumberKinds.o SphereGeom.o IntegerList.o OutputWriter.o Logger.o
PANEL_MPI_OBJS=oldStripack.o ssrfpack.o ParticlesEdgesPanels2.o TracerAndVorticityDistributions.o \
	RefineRemesh.o PanelVelocityParallel.o

%.o : %.f90
	$(FF) -c $(FF_FLAGS) $< `mpif90 -showme:compile`

%.o : %.f
	$(FF) -c $(FF_FLAGS) -real-size 64 $<

############################
#  Program executables	   #
############################
interpParticlesVTK2NCL.exe: InterpPassivePointsVTKToNCL.o $(BASE_OBJS)
	ifort -O2 -o interpParticlesVTK2NCL.exe InterpPassivePointsVTKToNCL.o $(BASE_OBJS)
findStratModelPotVort.exe: StratModelPotVort.o $(BASE_OBJS) TracerAndVorticityDistributions.o
	$(FF) $(FF_FLAGS) -o findStratModelPotVort.exe StratModelPotVort.o $(BASE_OBJS) TracerAndVorticityDistributions.o
solidBodyMPI.exe: solidBody.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o solidBodyMPI.exe solidBody.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
bveDeltaTConvMPI.exe: BVEDeltaTConv.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o bveDeltaTConvMPI.exe BVEDeltaTConv.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`	
tripolePanelsMPI.exe: TripolePanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o tripolePanelsMPI.exe TripolePanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
advectBlockMDirectMPI.exe: AdvectBlockMDirect.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectBlockMDirectMPI.exe advectBlockMDirect.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
advectBlockMMPI.exe: AdvectBlockM.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectBlockMMPI.exe AdvectBlockM.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
zonalMeanPanelsMPI.exe: ZonalMeanPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o zonalMeanPanelsMPI.exe ZonalMeanPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
stratospherePanelsMPI.exe: StratospherePanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o stratospherePanelsMPI.exe StratospherePanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
twoVortsPanelsMPI.exe: TwoVortsPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o twoVortsPanelsMPI.exe twoVortsPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`	
advectGaussHillsMPI.exe: NonDivAdvection.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectGaussHillsMPI.exe NonDivAdvection.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
advectGaussHillsDirectInterpMPI.exe: AdvectionDirectInterp.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectGaussHillsDirectInterpMPI.exe AdvectionDirectInterp.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
advectSlottedCDirectInterpMPI.exe: AdvectDirectSlottedC.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectSlottedCDirectInterpMPI.exe AdvectDirectSlottedC.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
advectSlottedCMPI.exe: AdvectSlottedC.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectSlottedCMPI.exe AdvectSlottedC.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
advectCosBellsMPI.exe: AdvectCosBells.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectCosBellsMPI.exe AdvectCosBells.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`	
advectCosBellsDirectMPI.exe: AdvectCosineBellsDirect.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o advectCosBellsDirectMPI.exe AdvectCosineBellsDirect.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
rh2PanelsMPI.exe: RH2Panels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o rh2PanelsMPI.exe RH2Panels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
manyVortsPanelsMPI.exe: ManyVortsPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o manyVortsPanelsMPI.exe ManyVortsPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
jetPanelsMPI.exe: JetPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o jetPanelsMPI.exe JetPanels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`
rh4DirectPanelsMPI.exe: DirectRH4Panels.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o rh4DirectPanelsMPI.exe DirectRH4Panels.o $(BASE_OBJS) $(PANEL_MPI_OBJS) \
	`mpif90 -showme:link`								
rh4PanelsMPI.exe: RHWavePanelsMPI.o $(BASE_OBJS) $(PANEL_MPI_OBJS) 
	$(FF) $(FF_FLAGS) -o $@ RHWavePanelsMPI.o $(BASE_OBJS) $(PANEL_MPI_OBJS) `mpif90 -showme:link`
gaussVortMPI.exe: GaussVort.o $(BASE_OBJS) $(PANEL_MPI_OBJS)
	$(FF) $(FF_FLAGS) -o $@ GaussVort.o $(BASE_OBJS) $(PANEL_MPI_OBJS) `mpif90 -showme:link`	

############################
#  Program object files	   #
############################
BVEDeltaTConv.o : BVEDeltaTConv.f90 $(BASE_OBJS)
InterpPassivePointsVTKToNCL.o: InterpPassivePointsVTKToNCL.f90 $(BASE_OBJS)
StratModelPotVort.o: StratModelPotVort.f90 $(BASE_OBJS) TracerAndVorticityDistributions.o
RHWavePanelsMPI.o: RHWavePanelsMPI.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
DirectRH4Panels.o: DirectRH4Panels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
RH2Panels.o : RH2Panels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
ManyVortsPanels.o: ManyVortsPanels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
JetPanels.o: JetPanels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
ZonalMeanPanels.o: ZonalMeanPanels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
TwoVortsPanels.o: TwoVortsPanels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
TripolePanels.o: TripolePanels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
StratospherePanels.o: StratospherePanels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
solidBody.o: solidBody.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
NonDivAdvection.o: NonDivAdvection.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
AdvectionDirectInterp.o: AdvectionDirectInterp.f90 	$(BASE_OBJS) $(PANEL_MPI_OBJS)
AdvectSlottedC.o: AdvectSlottedC.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
AdvectDirectSlottedC.o: AdvectDirectSlottedC.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
AdvectBlockM.o: AdvectBlockM.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
AdvectBlockMDirect.o: AdvectBlockMDirect.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
AdvectCosBells.o: AdvectCosBells.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
AdvectCosineBellsDirect.o: AdvectCosineBellsDirect.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
DirectRH4Panels.o: DirectRH4Panels.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)
GaussVort.o: GaussVort.f90 $(BASE_OBJS) $(PANEL_MPI_OBJS)

############################
# Modules 				   #
############################
NumberKinds.o : NumberKinds.f90
SphereGeom.o : SphereGeom.f90 NumberKinds.o
IntegerList.o : IntegerList.f90 NumberKinds.o
OutputWriter.o : OutputWriter.f90 NumberKinds.o
Logger.o : Logger.f90 OutputWriter.o NumberKinds.o
oldStripack.o : oldStripack.f
ssrfpack.o : ssrfpack.f
STRIPACK.o : STRIPACK.f90
VoronoiPanels.o : VoronoiPanels.f90 $(BASE_OBJS) STRIPACK.o ParticlesEdgesPanels2.o
ParticlesEdgesPanels2.o : ParticlesEdgesPanels2.f90 $(BASE_OBJS)
TracerAndVorticityDistributions.o : TracerAndVorticityDistributions.f90 ParticlesEdgesPanels2.o VoronoiPanels.o $(BASE_OBJS)
RefineRemesh.o : RefineRemesh.f90 TracerAndVorticityDistributions.o ParticlesEdgesPanels2.o $(BASE_OBJS)
PanelVelocityParallel.o : PanelVelocityParallel.f90 TracerAndVorticityDistributions.o ParticlesEdgesPanels2.o $(BASE_OBJS)

############################
# VTK executables		   #
############################
plotRelVortAndInitLat.exe: vtkPlotRelVortAndInitLat.o
	g++ -o plotRelVortAndInitLat.exe vtkPlotRelVortAndInitLat.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotRelVort.exe: vtkPlotRelVortFrames.o
	g++ -o plotRelVort.exe vtkPlotRelVortFrames.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotJetVoronoi.exe: vtkPlotJet2Panel.o
	g++ -o plotJetVoronoi.exe vtkPlotJet2Panel.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotJetAMR.exe: vtkPlotJet.o
	g++ -o plotJetAMR.exe vtkPlotJet.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotRH4Vort.exe: vtkPlotRelVortFrame2.o
	g++ -o plotRH4Vort.exe vtkPlotRelVortFrame2.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotRH4Voronoi.exe: vtkPlotRHWaveVoronoi.o 
	g++ -o plotRH4Voronoi.exe vtkPlotRHWaveVoronoi.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotMTracer.exe: vtkPlotMTracer.o
	g++ -o plotMTracer.exe vtkPlotMTracer.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotRH4LatTransport.exe: vtkPlotRH4LatTransport.o
	g++ -o plotRH4LatTransport.exe vtkPlotRH4LatTransport.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotMergVorts.exe: vtkPlotMergVorts.o
	g++ -o plotMergVorts.exe vtkPlotMergVorts.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotTwoVorts.exe: vtkPlotTwoVorts.o
	g++ -o plotTwoVorts.exe vtkPlotTwoVorts.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)
plotRH4IUTAM.exe: vtkIUTAMrh4.o
	g++ -o $@ vtkIUTAMrh4.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)	
plotGaussVortIUTAM.exe: vtkIUTAMGaussVort.o
	g++ -o $@ vtkIUTAMGaussVort.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)	
plotUnifGaussVortIUTAM.exe: vtkIUTAMGaussVortUnif.o
	g++ -o $@ vtkIUTAMGaussVortUnif.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)	
plotJetIUTAM.exe: vtkIUTAMJet.o
	g++ -o $@ vtkIUTAMJet.o $(VTK_INCLUDE) $(VTK_LIB_DIR) $(VTK_LIBS)	

############################
# VTK object files		   #
############################
vtkPlotRH4LatTransport.o: $< vtkPlotRH4LatTransport.cpp 
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)	
vtkPlotRelVortAndInitLat.o: vtkPlotRelVortAndInitLat.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkPlotRelVortFrames.o: vtkPlotRelVortFrames.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkPlotJet2Panel.o: vtkPlotJet2Panel.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkPlotJet.o: vtkPlotJet.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkPlotRelVortFrame2.o: vtkPlotRelVortFrame2.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkPlotMTracer.o: vtkPlotMTracer.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkPlotMergVorts.o: vtkPlotMergVorts.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkPlotTwoVorts.o: vtkPlotTwoVorts.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)		
vtkPlotRHWaveVoronoi.o: vtkPlotRHWaveVoronoi.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkIUTAMrh4.o: vtkIUTAMrh4.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)	
vtkIUTAMGaussVort.o: vtkIUTAMGaussVort.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)	
vtkIUTAMGaussVortUnif.o: vtkIUTAMGaussVortUnif.cpp	
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)
vtkIUTAMJet.o: vtkIUTAMJet.cpp
	g++ -c -Wno-deprecated -O2 $< $(VTK_INCLUDE)			