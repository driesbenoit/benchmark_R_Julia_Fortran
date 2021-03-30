! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! Draws from the posterior distribution of the Bayesian probit model
! with normal prior on the regression coefficients.

! Algorithm based on:
! Albert and Chib (1993). Bayesian Analysis of Binary and Polychotomous
! Response Data. Journal of the American Statistical Association, 
! 88(422), 669-679.

! Licence GLPv3 (or higher)

subroutine binprobbayes(y,x,b0,bb0,r,thin,n,k,betadraw)

implicit none

! Precision statement
integer, parameter :: dp = kind(1.0d0)

! Input arguments
integer, intent(in) :: r,thin,n,k
integer, dimension(n), intent(in) :: y 
real(dp), dimension(n,k), intent(in) :: x 
real(dp), dimension(k), intent(in) :: b0 
real(dp), dimension(k,k), intent(in) :: bb0 

! Output arguments
real(dp), dimension(r,k), intent(out) :: betadraw 

! Internal arguments
integer :: ok,i,ii
real(dp), dimension(k) :: b,beta
real(dp), dimension(k,k) :: bb
real(dp), dimension(n) :: ystar 

! Set starting values
beta = 0.0_dp
ystar = 0.0_dp

bb = bb0
call dpotrf('U',k,bb,k,ok)
call dpotri('U',k,bb,k,ok)
do i = 1,(k-1)
        do ii = (i+1),k
                bb(ii,i) = bb(i,ii)
        enddo
enddo

b = matmul(bb,b0)
bb = bb + matmul(transpose(x),x)

call dpotrf('U',k,bb,k,ok)
call dpotri('U',k,bb,k,ok)
do i = 1,(k-1)
        do ii = (i+1),k
                bb(ii,i) = bb(i,ii)
        enddo
enddo

! Start mcmc
do i = 1,r
        
        ! Draw new value for ystar
        do ii = 1,n
                if (y(ii) .eq. 0) then
                        call rtnorm(0.0_dp,.FALSE.,dot_product(x(ii,:),beta),1.0_dp,ystar(ii))
                else
                        call rtnorm(0.0_dp,.TRUE.,dot_product(x(ii,:),beta),1.0_dp,ystar(ii))
                endif
        enddo

        ! Draw new value for beta
        call rmvnorm(matmul(bb,(b+matmul(transpose(x),ystar))),bb,k,beta)

        ! Save draw
        betadraw(i,:) = beta

enddo

end subroutine binprobbayes 
