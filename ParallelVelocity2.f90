module ParallelVelocity2Module

use NumberKindsModule
use SphereGeomModule
use LoggerModule
use ParticlesEdgesPanelsModule
use TracerAndVorticityDistributionModule

implicit none

include 'mpif.h'

private
public totalKE
public kineticEnergy 	! KE at each panel
public BVERK4, AdvectionRK4, DivergentAdvectionRK4, StratosphereRK4 ! timestepping routines
public ResetRK4	 ! memory free routine
public ResetArea ! Set to true for divergent flow
public SetSmooth, GetSmooth
public InitializeMPIRK4, FinalizeMPIRK4, LoadBalance


end module
