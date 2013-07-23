module NumberKindsModule
	implicit none
	!
	! This module provides basic numerical data types and constants.
	!
	private
	public KREAL, KINT, KLOG, PI
		
	integer, parameter :: KREAL = kind(0.d0)
	integer, parameter :: KINT = kind(1)
	integer, parameter :: KLOG = kind(.TRUE.)
	real(KREAL), parameter :: PI = 3.141592653589793238462643383279502797479068098137295573004504331874296718662975536062731407582759857177734375d0
end module NumberKindsModule
