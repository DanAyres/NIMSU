program test

    use ORTHPOL_INTERFACE

    implicit none

    integer :: N_Coeff, SOLVER, N_Unions, Dist
    double precision :: Domain(2), Dist_args(2), scale
    
    double precision, allocatable   :: alpha(:), beta(:), PCoeffs(:,:)
    logical :: LeftInf, RightInf


    N_Coeff = 3
    SOLVER = 1
    N_Unions = 4
    Dist=1
    
    Domain(1)=0.0
    Domain(2) = 1.0
    scale = 1.0
    
    
    
    allocate(alpha(N_Coeff), beta(N_Coeff))

    call recursion_coefficients(N_Coeff, SOLVER, N_Unions, Domain,      &
                                LeftINF, RightINF, Dist, Dist_args,     &
                                scale, alpha, beta)

    allocate(PCoeffs(0:N_Coeff-1,0:N_Coeff-1) )
    call polynomial_coefficients(N_Coeff, alpha, beta, PCoeffs)

    print *, PCoeffs(0,:)
    print *, PCoeffs(1,:)
    print *, PCoeffs(2,:)

end program test
