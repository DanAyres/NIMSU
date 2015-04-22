module ProdDists

implicit none

contains

    double precision function Uniform(dx,i,Dist_args)
    !------------------------------------------------------------------------80
    !
    !   PDF of the uniform distribution on [-1:1]
    !
    !------------------------------------------------------------------------80
    double precision dx,Dist_args(2)
    integer :: i
!   Uniform = 0.5d0
    Uniform = 1.0d0
    return
    end function Uniform


    double precision function NormPDF(dx,i,Dist_args)
    !------------------------------------------------------------------------80
    !
    !   PDF of the normal distribution on [-oo:oo]
    !
    !------------------------------------------------------------------------80
    double precision dx,Dist_args(2)
    integer :: i
    NormPDF=dexp(-dx*dx)
    return
    end function NormPDF


    double precision function BetaPDF(dx,i,Dist_args)
    !------------------------------------------------------------------------80
    !
    !   PDF of the beta distribution on [-1:1]
    !
    !   Beta = (1-x)^alpha * (1+x)^beta
    !
    !------------------------------------------------------------------------80
    double precision dx,Dist_args(2)
    integer :: i
    BetaPDF = (1-dx)**Dist_args(1) * (1+dx)**Dist_args(2)
    return
    end function BetaPDF


    double precision function GamaPDF(dx,i,Dist_args)
    !------------------------------------------------------------------------80
    !
    !   PDF of the Gamma distribution on [0:oo]
    !
    !   Gama = x^(k-1) * exp ( x/theta ) 
    !
    !------------------------------------------------------------------------80
    double precision dx,Dist_args(2)
    integer :: i
    GamaPDF = dx**Dist_args(1) * dexp(-dx) 
    return
    end function GamaPDF

end module ProdDists
