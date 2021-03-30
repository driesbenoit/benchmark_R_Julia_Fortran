! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! This code generates one draw from the multivariate normal 
! distribution. 
! This code is based on the Cholesky decomposition of the
! variance-covariance matrix.

! Input arguments:
!	- mean   : mean vector
!	- varcov : varcov matrix 
!	- k      : dimension of mvnorm 

! Output arguments:
!	- fn_val	: random draw from N(0,1) distribution

! Licence: GPLv3 (or higher) 

subroutine rmvnorm(mean,varcov,k,fn_val)

implicit none

! Precision statement
integer, parameter :: dp = kind(1.0d0)

! Input arguments
integer, intent(in) :: k
real(dp), dimension(k), intent(in) :: mean
real(dp), dimension(k,k), intent(in) :: varcov

! Output arguments
real(dp), dimension(k), intent(out) :: fn_val 

! Internal arguments
integer :: ok,i1,i2
real(dp), dimension(k,k) :: cvarcov


! Cholesky root of varcov
cvarcov = varcov
call dpotrf('U',k,cvarcov,k,ok)
do i1 = 1,(k-1)
  do i2 = (i1+1),k
    cvarcov(i2,i1) = 0.0_dp
  enddo
enddo

! Standard normal vector
do i1 = 1,k
  call rnorm(fn_val(i1))
enddo

! Linear transformation
fn_val = mean + matmul(fn_val,cvarcov)

end subroutine rmvnorm
