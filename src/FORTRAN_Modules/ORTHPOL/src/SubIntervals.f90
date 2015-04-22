module SubInts

implicit none


contains


    subroutine Calculate_SubIntervals(N_Unions, Domain, LeftINF, RightINF,    &
                                      EndL, EndR, Scale)
    !------------------------------------------------------------------------80
    !
    !   Dan Ayres  11/2/2015
    !
    !   Given the domain and number of sub-intervals, calculate the 
    !   left and right edges.
    !    
    !
    !
    !   Args:
    !       N_Unions    : Number of sub-intervals                           int
    !       Domain      : Support of the polynomials                  double(2)
    !       LeftINF     : Left extreme is infinity                      logical
    !       RightINF    : Right extreme is infinity                     logical
    !       Scale       : Interval width for infinite domains           double
    !
    !   Returns:
    !       EndL        : Left boundaries                      double(N_Unions)
    !       EndR        : Right boundaries                     double(N_Unions)
    !
    !------------------------------------------------------------------------80
    
    integer             :: N_Unions
    double precision    :: Domain(2), Scale
    Logical             :: LeftINF, RightINF
    double precision    :: EndL(N_Unions), EndR(N_Unions)

    double precision    :: delta
    integer             :: i

    if (N_Unions==1) then
        EndL(1) = Domain(1)
        EndR(N_Unions) = Domain(2)
        return
    end if


    if (.not.LeftInf .and. .not.RightInf) then

        EndL(1) = Domain(1)
        EndR(N_Unions) = Domain(2)    
        delta = ( Domain(2) - Domain(1) ) / real(N_Unions)

    else if (LeftInf .and. .not.RightInf) then

        EndR(N_Unions) = Domain(2)    
        EndL(1) = Domain(2)-N_Unions*scale
        delta = scale

    else if (.not.LeftInf .and. RightInf) then

        EndR(N_Unions) = Domain(1) + N_Unions*scale    
        EndL(1) = Domain(1)
        delta = scale
        
    else if (LeftInf .and. RightInf) then
    
        EndL(1) = 0.0 - real(N_Unions)/real(2) * scale        
        EndR(N_Unions) = EndL(1) + N_Unions*scale    
        delta = scale
        
    end if


    do i=1,N_Unions-1
        EndL(i+1) =  EndL(i) + delta
        EndR(i) = EndL(i+1)
    end do            
    
    return
    end subroutine Calculate_SubIntervals

end module SubInts
