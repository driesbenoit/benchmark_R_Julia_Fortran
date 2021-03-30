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

program test_bayesprobit

        implicit none

        ! Precision statement
        integer, parameter :: dp = kind(1.0d0)

        integer :: i
        integer :: tic,toc,rate
        integer, dimension(5000) :: y
        !real(dp) :: tic, toc 
        real(dp), dimension(3) :: b0, bayesest 
        real(dp), dimension(3,3) :: bb0 
        real(dp), dimension(5000,3) :: x, betadraw

        ! Read in the data
        open(10, file="probitdata.csv",access='sequential',form="formatted")
        read(10,*) ! Skip first row

        do i = 1,5000
                read(10,*) y(i), x(i,:)
        enddo

        ! Prior
        b0 = 0.0_dp
        bb0 = 0.0_dp
        do i = 1,3
                bb0(i,i) = 1000.0_dp
        enddo

        ! Execute algorithm
        ! NOTE: impliciat conversion of y to logical
        !call cpu_time(tic)
        call system_clock(tic,rate)
        call binprobbayes(y,x,b0,bb0,5000,5000,3,betadraw)
        !call cpu_time(toc)
        call system_clock(count=toc)

        ! Write time
        !print *, 'Time in seconds:',toc-tic
        print *, 'Time in seconds:',(toc-tic)/real(rate)

        ! Write result 
        do i = 1,3
                bayesest(i) = sum(betadraw(:,i))/5000.0_dp
        enddo
        print *, 'Bayes estimate:', bayesest
end
