module mcdis_new

use quadrature
use Stieltjes

contains


    subroutine mcdis2(N_Coeff, ncapm, N_Unions, Solver, eps, beta, alpha, kount, ierr,&
                      ie, finl, finr, endl, endr, dwf, Dist_args)

    implicit none

!    Input:  N_Coeff    - - the number of recursion coefficients desired;
!                     type integer
!            ncapm  - a discretization parameter indicating an upper
!                     limit of the fineness of the discretization;
!                     ncapm=500  will usually be satisfactory; type
!                     integer
!            N_Unions  - -  the number of disjoint intervals in the
!                     continuous part of the spectrum; type integer    


!    Returns:  alpha,beta - arrays of dimension n, holding as k-th
!                     element  alpha(k-1), beta(k-1), k=1,2,...,n,
!                     respectively
!             ncap  - an integer indicating the fineness of the
!                     discretization that yields convergence within
!                     the eps-tolerance
!             kount - the number of iterations used
!             ierr  - an error flag, equal to  0  on normal return,
!                     equal to  -1  if  n  is not in the proper range,
!                     equal to  i  if there is an error condition in
!                     the discretization of the i-th interval,
!                     and equal to  ncapm  if the discretized 
!                     Stieltjes procedure does not converge within the
!                     discretization resolution specified by  ncapm
!             ie - -  an error flag inherited from the routine  sti
!                     or  lancz  (whichever is used)

    
    integer, intent(in)             :: N_Coeff, ncapm, N_Unions, Solver
    double precision, intent(in)    :: eps
    double precision, intent(inout) :: endl(N_Unions), endr(N_Unions)
    logical, intent(inout)          :: finl, finr
    
    integer, intent(out)    :: ierr, ie         ! Error
    integer, intent(out)    :: kount        ! Number of iterations used
    double precision, intent(out)   :: beta(N_Coeff), alpha(N_Coeff)
    double precision, external    :: dwf
    double precision, intent(inout) :: Dist_args(2)
    !    
    ! Working space
    !
    double precision :: be(N_Coeff)
    double precision :: xm(N_Unions*ncapm),wm(N_Unions*ncapm),&
                        p0(N_Unions*ncapm),p1(N_Unions*ncapm),&
                        p2(N_Unions*ncapm)
    double precision :: x(ncapm), w(ncapm), xfer(ncapm), wfer(ncapm)
    !
    ! Locals    
    !
    integer :: incap, ncap
    integer :: k, i, iterator
    integer :: mtncap, im1tn
    logical :: flag
    integer, parameter :: Max_Its = 30
    
    
    
        
    if(N_Coeff.lt.1) then
        ierr=-1
        return
    end if

!
!   Initialization
!            
    incap=1
    kount=-1
    ierr=0
    
    do k=1,N_Coeff
        beta(k)=0.d0
    end do
    ncap=(2*N_Coeff-1)
    
    do iterator=1,Max_Its
    
        do k=1,N_Coeff
            be(k)=beta(k)
        end do
        kount=kount+1
        if(kount.gt.1) incap=2**(kount/5)*N_Coeff
        ncap=ncap+incap
        if(ncap.gt.ncapm) then
            ierr=ncapm
            return
        end if
        !
        !   Discretization of the inner product
        !
        mtncap=N_Unions*ncap
        do i=1,N_Unions
            im1tn=(i-1)*ncap
            call dqgp(ncap,x,w,i,ierr,N_Unions,finl,finr,endl,endr,xfer,wfer,dwf, Dist_args)
            
            if(ierr.ne.0) then
              ierr=i
              return
            end if
            
            do k=1,ncap
              xm(im1tn+k)=x(k)
              wm(im1tn+k)=w(k)
            end do
        end do
        !
        !   Computation of the desired recursion coefficients
        !   
        if(Solver.eq.1) then
            call dsti(N_Coeff,mtncap,xm,wm,alpha,beta,ie,p0,p1,p2)
        else
            call dlancz(N_Coeff,mtncap,xm,wm,alpha,beta,ie,p0,p1)
        end if
        !
        !   Convergence check
        !
        flag=.true.
        do k=1,N_Coeff
            if(dabs(beta(k)-be(k)).gt.eps*dabs(beta(k))) flag=.false.
        end do
        
        if (flag) exit
    
    end do
    
    return

    end subroutine mcdis2
    
    
    
end module mcdis_new
