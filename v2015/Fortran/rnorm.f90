! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! This code generates one draw from the standard normal 
! distribution. Note that more efficient code is possible
! when more than one normal draw is required.
! This code is based on the Box-Muller method.

! Output arguments:
!	- fn_val	: random draw from N(0,1) distribution

! Licence: GPLv3 (or higher) 


subroutine rnorm(fn_val)

implicit none

! precision statement:
integer, parameter :: dp = kind(1.0d0)

! output arguments:
real(dp), intent(out) :: fn_val

! internal arguments:
real(dp) :: pi
real(dp), dimension(1:2) :: u

pi = 3.14159265358979323846_dp

call random_number(u)

fn_val = sqrt(-2.0_dp*log(u(1))) * cos(2.0_dp*pi*u(2))

end subroutine rnorm
