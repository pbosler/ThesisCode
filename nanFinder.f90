program nanFinder

use NumberKindsModule
use SphereGeomModule
use TracerAndVorticityDistributionModule

implicit none

double precision :: one, zero

one = 1.0d0
zero = 0.0d0


print '(A,F24.18)', "ArcTan(1/0) = ",atan2(one,zero)

print '(A,F24.18)', "A(0) = ", Juckes_A(zero)
print '(A,F24.18)', "A(1) = ", Juckes_A(one)

print '(A,F24.18)', "B(0) = ", Juckes_B(zero)

print '(A,F24.18)', "Longitude( [1,0,0]) = ", Longitude([one,zero,zero])
print '(A,F24.18)', "Latitude([1,0,0]) = ",Latitude([one,zero,zero])

print '(A,F24.18)', "JuckesForcing([0,0,1],0) = ", Juckes_Forcing([zero,zero,one],zero)
print '(A,F24.18)', "JuckesForcing([0,0,1],1) = ", Juckes_Forcing([zero,zero,one],one)
print '(A,F24.18)', "JuckesForcing([1,0,0],0) = ", Juckes_Forcing([one,zero,zero],zero)
print '(A,F24.18)', "JuckesForcing([1,0,0],1) = ", Juckes_Forcing([one,zero,zero],one)

end program
