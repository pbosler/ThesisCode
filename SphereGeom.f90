! Pete Bosler
! October 26, 2011
! This module defines functions for performing geometric calculation 
! on the surface of the unit sphere.
! filename : SphereGeom.f90
! Last modified : 11-9-12
module SphereGeomModule

use NumberKindsModule

implicit none

public

interface SphereTriArea
	module procedure SphereTriAreaVector
	module procedure SphereTriAreaComponents
end interface

interface SphereDistance
	module procedure SphereDistanceVector
	module procedure SphereDistanceComponents
end interface

interface Latitude
	module procedure LatitudeVector
	module procedure LatitudeComponents
end interface

interface Longitude
	module procedure LongitudeVector
	module procedure LongitudeComponents
	module procedure LongitudeComponents2
end interface

contains

!----------------
! Basic geometry : length, area, coordinates, etc.
!----------------

function ChordDistance(xyzA, xyzB)
	! Outputs the Euclidean distance between two points in R3.
	real(kreal), intent(in) :: xyzA(3), xyzB(3)
	real(kreal) :: ChordDistance
	ChordDistance = sqrt( sum( (xyzB - xyzA)*(xyzB-xyzA) ) )
end function


function SphereDistanceVector( xyzA, xyzB )
! Finds the great circle distance between two points (xyzA and xyzB) on the unit sphere
! Uses atan formulation (rather than acos) for numeric stability for nearby input points.
	! Calling parameters
	real(KREAL), dimension(3) :: xyzA, xyzB
	real(KREAL) :: SphereDistanceVector
	! Local variables
	real(KREAL), dimension(3) :: crossProd
	real(KREAL) :: dotProd , crossNorm
	
	crossProd = [xyzA(2)*xyzB(3)-xyzB(2)*xyzA(3),xyzB(1)*xyzA(3)-xyzA(1)*xyzB(3),&
    xyzA(1)*xyzB(2)-xyzB(1)*xyzA(2) ]
    
    crossNorm = dsqrt(sum(crossProd*crossProd))
    
    dotProd = xyzA(1)*xyzB(1)+xyzA(2)*xyzB(2)+xyzA(3)*xyzB(3)
    
    SphereDistanceVector = datan2(crossNorm,dotProd)
end function


function SphereDistanceComponents(xA, yA, zA, xB, yB, zB)
	real(kreal), intent(in) :: xA, yA, zA
	real(kreal), intent(in) :: xB, yB, zB
	real(kreal) :: SphereDistanceComponents
	real(kreal) :: cp1, cp2, cp3, cpNorm, dp
	
	cp1 = yA*zB - yB*zA
	cp2 = xB*zA - xA*zB
	cp3 = xA*yB - xB*yA
	
	cpNorm = sqrt( cp1*cp1 + cp2*cp2 + cp3*cp3)
	
	dp = xA*xB + yA*yB + zA*zB
	
	SphereDistanceComponents = atan2(cpNorm,dp)
end function




function SphereMidpoint(xyzA, xyzB)
! Finds the midpoint of two points on the unit sphere by finding the midpoint of the chord
! connecting the two points, then projecting the chord midpoint to the sphere.
	! Calling parameters
	real(KREAL), dimension(3) :: xyzA, xyzB
	real(KREAL), dimension(3) :: SphereMidpoint
	
	SphereMidpoint = (xyzA + xyzB)/2.0_KREAL
	SphereMidpoint = SphereMidpoint/sqrt(sum(sphereMidpoint*sphereMidpoint))
end function 

function SphereTriCenter(xyzA, xyzB, xyzC)
! Finds the midpoint of three points on the unit sphere by find their average position in Cartesian
! coordinates, then projecting that average onto the sphere.
	! Calling parameters
	real(KREAL), dimension(3) :: xyzA, xyzB, xyzC
	real(KREAL), dimension(3) :: SphereTriCenter

	SphereTriCenter = (xyzA + xyzB + xyzC)/3.0_KREAL
	SphereTriCenter = SphereTriCenter/sqrt(sum(sphereTriCenter*sphereTriCenter))
end function


function SphereQuadCenter(xyzA, xyzB, xyzC, xyzD)
! Finds the midpoint of four points on the unit sphere by finding their average position in 
! Cartesian coordinates, then projecting that average onto the sphere.
	! Calling parameters
	real(KREAL), dimension(3) :: xyzA, xyzB, xyzC, xyzD
	real(KREAL), dimension(3) :: sphereQuadCenter
	
	SphereQuadCenter = (xyzA + xyzB + xyzC + xyzD)/4.0_KREAL
	SphereQuadCenter = SphereQuadCenter/sqrt(sum(sphereQuadCenter*sphereQuadCenter))
end function


function SphereTriAreaVector(xyzA, xyzB, xyzC)
! Calculates the area of a spherical triangle on the unit sphere
!   NOTE : This function requires function sphereDistance.
	! Calling parameters
	real(KREAL), dimension(3) :: xyzA, xyzB, xyzC
	real(KREAL) :: SphereTriAreaVector
	! Local variables
	real(KREAL) :: side1, side2, side3, halfPerimeter, zz
	
	side1 = SphereDistance(xyzA,xyzB)
	side2 = SphereDistance(xyzB,xyzC)
	side3 = SphereDistance(xyzC,xyzA)
	
	halfPerimeter = (side1 + side2 + side3)/2.0_KREAL
	
	zz = dtan(halfPerimeter/2.0_KREAL)*dtan( (halfPerimeter-side1)/2.0_KREAL )*&
		dtan( (halfPerimeter - side2)/2.0_KREAL )*dtan( (halfPerimeter - side3)/2.0_KREAL )
	
	SphereTriAreaVector = 4.0_KREAL * datan2(sqrt(zz),1.0_KREAL)
end function


function SphereTriAreaComponents(xa,ya,za, xb,yb,zb, xc,yc,zc)
	real(kreal), intent(in) :: xa,ya,za
	real(kreal), intent(in) :: xb,yb,zb
	real(kreal), intent(in) :: xc,yc,zc
	real(kreal) :: SphereTriAreaComponents
	real(kreal) :: s1,s2,s3, halfPerim, zz
	
	s1 = SphereDistanceComponents(xa,ya,za,xb,yb,zb)
	s2 = SphereDistanceComponents(xb,yb,zb,xc,yc,zc)
	s3 = SphereDistanceComponents(xc,yc,zc,xa,ya,za)
	
	halfPerim = (s1+s2+s3)/2.0_kreal
	zz = tan( halfPerim/2.0_kreal)*tan( (halfPerim-s1)/2.0_kreal) * &
		 tan( (halfPerim-s2)/2.0_kreal)*tan( (halfPerim-s3)/2.0_kreal)
	
	SphereTriAreaComponents = 4.0_kreal*atan(sqrt(zz))
end function



function PlaneTriArea(xyzA,xyzB,xyzC)
	! Outputs the area of a planar triangle in R3 with vertices xyzA,B,C.
	real(kreal) :: PlaneTriArea
	real(kreal), intent(in) :: xyzA(3), xyzB(3), xyzC(3)
	PlaneTriArea = crossMagnitude( xyzB-xyzA, xyzC-xyzA)/2.0_kreal
end function


function longitudeVector(xyz)	
	! Outputs the longitude of a point on the unit sphere.
	real(KREAL) :: longitudeVector
	real(KREAL), intent(in) :: xyz(3)
	longitudeVector = atan4(xyz(2),xyz(1))
end function


function LongitudeComponents(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: LongitudeComponents
	LongitudeComponents = atan4(y,x)
end function


function LongitudeComponents2(x,y)
	real(kreal), intent(in) :: x, y
	real(kreal) :: LongitudeComponents2
	LongitudeCOmponents2 = atan4(y,x)
end function


function latitudeVector(xyz)
	! Outputs the latitude of a point on the unit sphere.
	real(KREAL) :: latitudeVector
	real(KREAL), intent(in) :: xyz(3)
	latitudeVector = datan2(xyz(3),dsqrt(xyz(1)**2 + xyz(2)**2))
end function


function LatitudeComponents(x,y,z)
	real(kreal), intent(in) :: x, y, z
	real(kreal) :: LatitudeComponents
	LatitudeComponents = atan2(z,sqrt(x*x+y*y))
end function


function atan4(y,x)
	!This function computes the inverse tangent (like atan2) but outputs angles in the range
	! 0 to 2 pi (rather than -pi to pi).
	! Adapted from John Burkhardt: http://people.sc.fsu.edu/~jburkhardt/m_src/halton/atan4.m
	real(kreal), intent(in) :: y,x
	real(kreal) :: atan4
	real(kreal) :: absY, absX, theta
	if ( x == 0.0_kreal) then
		if ( y > 0.0_kreal) then
			atan4 = PI/2.0_kreal
		elseif ( y < 0.0_kreal) then
			atan4 = 3.0_kreal*PI/2.0_kreal
		elseif ( y == 0.0_kreal) then
			atan4 = 0.0_kreal
		endif
	elseif ( y == 0.0_kreal) then
		if ( x > 0.0_kreal) then
			atan4 = 0.0_kreal
		elseif ( x < 0.0_kreal) then
			atan4 = PI
		endif
	else
		absY = abs(y)
		absX = abs(x)
		theta = datan2(absY,absX)
		if ( (x>0.0_kreal) .and. (y>0.0_kreal)) then
			atan4 = theta
		elseif ( (x < 0.0_kreal) .and. (y > 0.0_kreal)) then
			atan4 = pi -theta
		elseif ( (x < 0.0_kreal) .and. (y < 0.0_kreal)) then
			atan4 = pi+ theta
		elseif ( (x > 0.0_kreal) .and. (y<0.0_kreal)) then
			atan4 = 2.0_kreal*PI - theta
		endif		
	endif		
end function


function crossMagnitude(xyzA,xyzB)
	! Computes the magnitude of xyzA cross xyzB
	real(kreal) :: crossMagnitude
	real(kreal), intent(in) :: xyzA(3), xyzB(3)
	crossMagnitude = sqrt( (xyzA(2)*xyzB(3) - xyzA(3)*xyzB(2))*(xyzA(2)*xyzB(3) - xyzB(2)*xyzA(3)) + &
	   					   (xyzA(3)*xyzB(1) - xyzA(1)*xyzB(3))*(xyzA(3)*xyzB(1) - xyzA(1)*xyzB(3)) + &
	   					   (xyzA(1)*xyzB(2) - xyzA(2)*xyzB(1))*(xyzA(1)*xyzB(2) - xyzA(2)*xyzB(1)))
end function


function crossProduct(xyzA,xyzB)
	! Computes the cross product vector xyzA cross xyzB
	real(kreal) :: crossProduct(3)
	real(kreal), intent(in) :: xyzA(3), xyzB(3)
	crossProduct(1) = xyzA(2)*xyzB(3) - xyzA(3)*xyzB(2)
	crossProduct(2) = xyzA(3)*xyzB(1) - xyzA(1)*xyzB(3)
	crossProduct(3) = xyzA(1)*xyzB(2) - xyzA(2)*xyzB(1)
end function	

function Determinant(xA,xB,xC)
	! Computes the vector determinant (triple product) of xA, xB, and xC
	! This result will be positive if xA lies to the left of the directed arc
	! xB to xC, and negative if xA lies to the right of the arc xB->xC.  
	real(kreal) :: Determinant
	real(kreal), intent(in) :: xA(3), xB(3), xC(3)
	real(kreal) :: cross(3)
	cross = CrossProduct(xB,xC)
	Determinant = sum(xA*cross)
end function


function RandomSpherePoint()
	! Outputs a random point on the unit sphere.
	! Note : the seed is not initialized, so each sequence of 
	!        random points will be the same.
	real(kreal) :: RandomSpherePoint(3)
	real(kreal) :: x, y, z
	!call Random_Seed()
	call Random_Number(x)
	call Random_Number(y) 
	call Random_Number(z)
	x = -1.0_kreal + 2.0_kreal*x
	y = -1.0_kreal + 2.0_kreal*y
	z = -1.0_kreal + 2.0_kreal*z
	RandomSpherePoint(1) = x
	RandomSpherePoint(2) = y
	RandomSpherePoint(3) = z
	RandomSpherePoint = RandomSpherePoint/sqrt(sum(RandomSpherePoint*RandomSpherePoint))
end function

end module SphereGeomModule
