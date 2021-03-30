! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! Returns one draw from the truncated normal distribution

! Algorithm based on:
! Geweke, J. (1991). Efficient Simulation From the Multivariate Normal 
! and Student t-Distributions Subject to Linear Constraints, in Computer 
! Sciences and Statistics Proceedings of the 23d Symposium on the 
! Interface, pp. 571-578.

! This subroutine makes use of the subroutines:
!	- rnorm		: Box-Muller method for random normal draws

! Input arguments:
! a             -	trucation point
! lb	        -	logical:        if .TRUE. then trucation (a,+Inf)
!		                    		if .FALSE. then truncation (-Inf,a)
! mu	        -	mean of trunc normal
! sigma         -	sd of trunc normal
! fn_val        -	random draw from trunc normal

! Licence: GPLv3 (or higher)


subroutine rtnorm(a, lb, mu, sigma, fn_val)
  
implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
logical, intent(in) :: lb
real(dp), intent(in) :: a, mu, sigma

! Output arguments:
real(dp), intent(out) :: fn_val

! Internal arguments:
real(dp) :: z, phi_z, az, c
real(dp), dimension(1:2) :: u


! Rescale truncation point
az=(a-mu)/sigma

if (lb) then
    c=az
  else 
    c=-az
endif

if (c<.45_dp) then

  ! normal rejection sampling
  do
    call rnorm(u(1))
    
    if (u(1)>c) exit
  end do
  z=u(1)

else

  ! exponential rejection sampling
  do
    ! Create exponential random variate z
    ! from uniform random variate u(1)
    call random_number(u)
    z = -log(u(1))/c

    phi_z = exp(-.5_dp * z**2_dp) !see Geweke
    if (u(2)<phi_z) exit
  end do
  z=z+c

end if

if (lb) then 
  fn_val = mu + sigma*z
else
  fn_val = mu - sigma*z
end if

end subroutine rtnorm
