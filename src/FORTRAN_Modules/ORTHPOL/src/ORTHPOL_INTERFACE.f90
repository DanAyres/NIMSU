

module ORTHPOL_INTERFACE
use mcdis_new
use ProdDists
use SubInts

implicit none

contains


    subroutine Normalisation_Factors(N_Poly, N_Quad, PCoeffs, abcissa,        &
                                     weights, NFactors)
    !------------------------------------------------------------------------80
    !
    !   Compute the normalisation factors of the polynomials.
    !
    !       <P_n(x)^2>
    !
    !   Args:
    !       N_Poly     : Number of polynomials                              int
    !       N_Quad     : Number of quadrature points                        int
    !       PCoeffs    : Coefficients of x powers     double(0:N_Poly,0:N_Poly)
    !       abcissa    :                                         double(N_Quad)
    !       weights    :                                         double(N_Quad)
    !
    !   Returns:
    !       NFactors   : Normalisation factors indexed from 0  double(o:N_Poly)
    !
    !------------------------------------------------------------------------80
    integer           :: N_Poly, N_Quad
    double precision  :: PCoeffs(0:N_Poly, 0:N_Poly) 
    double precision  :: abcissa(N_Quad), weights(N_Quad)
    
    double precision  :: NFactors(0:N_Poly)
    
    integer           :: ip, iq
    double precision  :: x(0:N_Poly), p(0:N_Poly, N_Quad)

!f2py intent(in) N_Poly, PCoeffs, abcissa, weights
!f2py intent(out) NFactors

    p=0.0d0    

    do iq=1, N_Quad

        ! Compute x^i    
        do ip=0, N_Poly
            x(ip) = abcissa(iq)**ip
        end do
        
        ! Compute P(x)
        do ip=0, N_Poly
            p(ip,iq) = sum(PCoeffs(:,ip)*x)
        end do
        
    end do
    
    ! Integrate
    NFactors = 0.d0
    do iq=1, N_Quad
        NFactors = NFactors + weights(iq) * p(:,iq)**2
    end do
    
    end subroutine Normalisation_Factors



!    double precision  function legmonic(n)
!    ! Normalisation factor for the monic Legendre polynomials
!        integer :: n
!        legmonic = (2.d0/(2.d0*n + 1.d0))*( (2**n * fact1(n)**2)/fact1(2*n) )**2
!    end

!    double precision function fact1(n)
!        integer :: n
!        double precision vals(0:10)
!        vals = (/1.d0, 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0, 40320.d0, 362880.d0, 3628800.d0/)
!        fact1 = vals(n)    
!    end 

    subroutine polynomial_coefficients(N_Coeff, alpha, beta, PCoeffs)
    !------------------------------------------------------------------------80
    !
    !   Calculate the coefficients of the powers of x in Pn(x) from
    !   the recursion coefficients.
    !
    !   Args:
    !       N_Coeff     : Number of recursion coefficients                  int
    !       alpha       : alpha recursion coefficients          double(N_Coeff)
    !       beta        : beta recursion coefficients           double(N_Coeff)
    !
    !   Returns:
    !       PCoeffs     : Coefficients of x powers  double(N_Coeff+1,N_Coeff+1)
    !
    !
    !------------------------------------------------------------------------80
    
    integer           :: N_Coeff
    double precision  :: alpha(N_Coeff), beta(N_Coeff)
    
    
    double precision  :: PCoeffs(0:N_Coeff-1, 0:N_Coeff-1) !(Xord,Pord)

!f2py intent(in) N_Coeff, alpha, beta
!f2py intent(out) PCoeffs
    
    integer :: Pord
    
    PCoeffs=0.d0
    PCoeffs(0,0) = 1.d0
    PCoeffs(1,1) = PCoeffs(0,0)
    PCoeffs(0,1) = -alpha(1)*PCoeffs(0,0)
    do Pord=2,N_Coeff-1
        PCoeffs(1:Pord,Pord) = PCoeffs(0:Pord-1,Pord-1)     ! x * Pn
        PCoeffs(0:Pord-1,Pord) = PCoeffs(0:Pord-1,Pord) &
                 - alpha(Pord) * PCoeffs(0:Pord-1,Pord-1)   ! -alpha_k * Pn
        PCoeffs(0:Pord-2,Pord) = PCoeffs(0:Pord-2,Pord) &
                 - beta(Pord) * PCoeffs(0:Pord-2,Pord-2)    ! -beta_k * Pn-1
    end do


    
    end subroutine polynomial_coefficients
    
    

    subroutine recursion_coefficients(N_Coeff, SOLVER, N_Unions, Domain,      &
                                      LeftINF, RightINF, Dist, Dist_args,     &
                                      scale, alpha, beta)
    !------------------------------------------------------------------------80
    !
    !   Dan Ayres  11/2/2015
    !
    !
    !   Computes the recursion coefficients, a_k and b_k, for an arbitrary 
    !   monic polynomial. The recursion relation is written as
    !
    !   P(k+1)(x) = (x-a(k))P(k)(x) - b(k)P(k-1)
    !
    !
    !   Two error flags  ierr, ie  are used which signal the occurrence 
    !   of an error condition in the quadrature process, or in the routine 
    !   sti  or  lancz  (whichever is used), respectively.
    !
    !
    !   Args:
    !       N_COEFF     : Number of recursion coefficients                  int
    !       SOLVER      : Soluition algorithm. 1=Stieltjes                 
    !                                          2=Lanczos.                   int
    !       N_UNIONS    : Number of sub-intervals                           int
    !       Domain      : Support of the polynomials                  double(2)
    !       LeftINF     : Left extreme is infinity                      logical
    !       RightINF    : Right extreme is infinity                     logical
    !       Dist        : Type of probability distribution                  int
    !                     1= Uniform pdf=1.0
    !                     2= Normal  pdf= exp(-x*x)
    !                     3= Beta pdf=(1-x)**Dist_args(1) * (1+x)**Dist_args(2)
    !                     4= Gamma
    !       Dist_args   : Parameter for the PDF                       double(2)
    !       scale       : interval size for infinite domain              double
    !
    !
    !   Returns:
    !       alpha       : alpha recursion coefficients          double(N_Coeff)
    !       beta        : beta recursion coefficients           double(N_Coeff)
    !    
    !
    !------------------------------------------------------------------------80
    implicit none
    
    
    integer            :: N_Coeff, SOLVER, N_Unions, Dist
    double precision   :: Domain(2), scale
    double precision   :: Dist_args(2)
    logical            :: LeftINF, RightINF

    double precision   :: alpha(N_Coeff), beta(N_Coeff)
    
!f2py intent(in) N_Coeff, SOLVER, N_Unions, Domain, Dist_args, LeftINF, RightINF
!f2py intent(out) alpha, beta
    
    double precision   :: EndL(N_Unions), EndR(N_Unions), eps
    
    
    integer     :: ncapm, kount, ierr, ie
    
    


    call Calculate_SubIntervals(N_Unions, Domain, LeftINF, RightINF, EndL,    &
                                EndR, Scale)    
    
    LeftINF = .not.LeftINF
    RightINF = .not.RightINF
    ncapm = 500
    eps=1000.d0*epsilon(eps)
    
    select case(Dist)
    
    case(1) ! Uniform
        call mcdis2(N_Coeff, ncapm, N_Unions, Solver, eps, beta, alpha, kount,&
                  ierr, ie, LeftINF, RightINF, EndL, EndR, Uniform, Dist_args)       

    case(2) ! Normal
        call mcdis2(N_Coeff, ncapm, N_Unions, Solver, eps, beta, alpha, kount,&
                  ierr, ie, LeftINF, RightINF, EndL, EndR, NormPDF, Dist_args)

    case(3) ! Beta
        call mcdis2(N_Coeff, ncapm, N_Unions, Solver, eps, beta, alpha, kount,&
                  ierr, ie, LeftINF, RightINF, EndL, EndR, BetaPDF, Dist_args)

    case(4) ! Gamma
        call mcdis2(N_Coeff, ncapm, N_Unions, Solver, eps, beta, alpha, kount,&
                  ierr, ie, LeftINF, RightINF, EndL, EndR, GamaPDF, Dist_args)
    case default

        stop
    end select
    
    end subroutine recursion_coefficients
    
    
    
    subroutine quadrature_rule(N_Coeff, alpha, beta, eps, abcissa, weights)
    !------------------------------------------------------------------------80
    !
    !   Calculate the Gauss quadrature scheme from the polynomial
    !   recursion coefficients.
    !
    !   Args:
    !       N_Coeff     : Number of recursion coefficients                  int
    !       alpha       : alpha recursion coefficients          double(N_Coeff)
    !       beta        : beta recursion coefficients           double(N_Coeff)
    !       eps         : the relative accuracy desired 
    !                     in the abcissa and weights                     double
    !
    !   Returns:
    !       abcissa     :                                       double(N_Coeff)
    !       weights     :                                       double(N_Coeff)
    !
    !
    !------------------------------------------------------------------------80    

    integer           :: N_Coeff
    double precision  :: alpha(N_Coeff), beta(N_Coeff), eps
    
    double precision  :: abcissa(N_Coeff), weights(N_Coeff)
    
    integer           :: ierr
    double precision  :: work(N_Coeff)
    
!f2py intent(in) N_Coeff, alpha, beta, eps
!f2py intent(out) abcissa, weights

    call dgauss(N_Coeff, alpha, beta, eps, abcissa, weights, ierr, work)
    
    end subroutine quadrature_rule


end module ORTHPOL_INTERFACE
