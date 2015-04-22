module GerstnerGriebel

implicit none


contains


    subroutine Init_Adaptive_Grid(NumDim, Active, Old, Indicator, N_Forward, N_Backward, Idx)
    !
    !   Initialise the adaptive sparse grid. The maximum level is equal to one.
    !
    !   3D example:
    !
    !   I   Idx       N_For     N_Back    Indicator
    !   1   0  0  0   2  3  4   -1 -1 -1   0
    !   2   1  0  0  -1 -1 -1    1 -1 -1   1
    !   3   0  1  0  -1 -1 -1   -1  1 -1   1
    !   4   0  0  1  -1 -1 -1   -1 -1  1   1
    !
    !   Args:
    !       NumDim
    !   Return:
    !       Active      : The active index locations
    !       Old         : The old index locations
    !       Indicator   : ==1:Active, ==0, Old
    !       N_Forward   : Forward neighbours
    !       N_Backward  : Backward neighbours
    !       Idx         : The actual indexes
    !
    implicit none
    
    integer     :: NumDim, MaxVals, Active(NumDim), Old(1), Indicator(NumDim+1)
    integer     :: N_Forward(NumDim,NumDim+1), N_Backward(NumDim,NumDim+1), Idx(NumDim,NumDim+1)
    
    integer     :: i
    
    !f2py intent(in)  NumDim, MaxVals
    !f2py intent(out) Active, Old, Indicator, N_Forward, N_Backward, Idx
    
    N_Forward = -1
    N_Backward=-1
    Idx = 0
    Old(1) = 1
    Indicator(1)=0

    do i=1, NumDim
        Active(i) = i+1
        Indicator(i+1) = 1
        N_Forward(i,1) = i+1
        N_Backward(i,i+1) = 1
        Idx(i,i+1) = 1
    end do
    
    end subroutine Init_Adaptive_Grid



    subroutine CalculateNeighbours(NumDim, N_Backward, N_Forward, Active_Max, N_Idx, Idx, Indicator, Level_Max)
    !
    !   Calculate all of the admissible forward neighbours for an index set at position Active_Max.
    !
    !   Args:
    !       NumDim
    !       Active_Max
    !   Return:
    !       N_Forward   : Forward neighbours
    !       N_Backward  : Backward neighbours
    !       N_Idx       : Total number of indexes
    !       Idx         : The actual indexes
    !       Indicator   : ==1:Active, ==0, Old
    !
    implicit none
    
    integer     :: NumDim, Active_Max, N_Idx, Level_Max
    integer     :: N_Backward(:,:), N_Forward(:,:), Idx(:,:), Indicator(:)

    !f2py intent(in)  NumDim, Active_Max, Level_Max
    !f2py intent(in,out) N_Backward, N_Forward, Idx, Indicator, N_Idx
    
    integer     :: kb, kbf
    integer     :: dim1, dim2, dim3
    integer     :: mloc(1)
    
    
    ! If all other dimensions are zero, the grid is allowed to grow up the axis
    if ( sum(Idx(:,Active_max)) == maxval(Idx(:,Active_max)) ) then
    
        mloc= maxloc(Idx(:,Active_max))
        dim1=mloc(1)
        
        if ( Idx(dim1,Active_max) < Level_Max ) then
        
            N_Idx = N_Idx + 1
            Idx(:,N_Idx) = Idx(:,Active_max)
            Idx(dim1,N_Idx) = Idx(dim1,N_Idx) + 1

            ! Neighbours
            N_Forward(:,N_Idx) = -1
            N_Backward(:,N_Idx) = -1
            N_Backward(dim1,N_Idx) = Active_max
            N_Forward(dim1,Active_max) = N_Idx
            
            Indicator(N_Idx) = 1

        end if        
        
    end if
    
    do dim1=1, NumDim
        
        
        !   Example:
        !    ______
        !   |      |      
        !   |  k -> i n  
        !   |__^___|__^___
        !   |  j   |  j   |
        !   |  kb -> i kbf ?
        !   |______|______|
        dim2 = i4_wrap ( dim1 + 1, 1, NumDim ) ! Calculate direction j
        

        if (Idx(dim2,Active_max)==0) cycle

        kb  = N_Backward(dim2,Active_max)         ! Backward neighbour in direction j
        kbf = N_Forward(dim1,kb)               ! Forward neighbour in direction i
            
        ! n is not a valid forward neighbour of k in direction i
        if ( kbf == -1 ) then
            cycle
        end if
        ! If kbf has a forward neighbour then the indices already exist
        if( N_Forward(dim2,kbf) .ne. -1 ) then
            cycle
        end if
        ! kbf cannot be in the active set
        if (Indicator(kbf) == 1) then
            cycle
        end if
        
        ! No dimension of the new index can be greater than Level_Max
        
        if ( Idx(dim1,Active_max) .ge. Level_Max ) then
            cycle
        end if
        
        ! Add the new indices
        N_Idx = N_Idx + 1
        Idx(:,N_Idx) = Idx(:,Active_max)
        Idx(dim1,N_Idx) = Idx(dim1,N_Idx) + 1
        
        ! Neighbours
        N_Forward(:,N_Idx) = -1
        
        N_Backward(:,N_Idx) = -1
        N_Backward(dim1,N_Idx) = Active_max
        N_Forward(dim1,Active_max) = N_Idx
        
        N_Backward(dim2,N_Idx) = kbf
        N_Forward(dim2,kbf)  = N_Idx
        
        ! Check all other dimensions for Neighbours
        do dim2 = dim1 + 2, dim1 + NumDim - 1

            dim3 = i4_wrap ( dim2, 1, NumDim )
            if (Idx(dim3,Active_max)==0) then
                cycle
            end if
            kb  = N_Backward(dim3,Active_max)
            kbf = N_Forward(dim1,kb)

            if ( kbf == -1 ) then
                cycle
            end if

            N_Backward(dim3,N_Idx) = kbf
            N_Forward(dim3,kbf) = N_Idx

        end do
        
        Indicator(N_Idx) = 1
        
     end do
     
     
    end subroutine CalculateNeighbours


    function i4_modp ( i, j )

    !*****************************************************************************80
    !
    !! I4_MODP returns the nonnegative remainder of I4 division.
    !
    !  Discussion:
    !
    !    If
    !      NREM = I4_MODP ( I, J )
    !      NMULT = ( I - NREM ) / J
    !    then
    !      I = J * NMULT + NREM
    !    where NREM is always nonnegative.
    !
    !    The MOD function computes a result with the same sign as the
    !    quantity being divided.  Thus, suppose you had an angle A,
    !    and you wanted to ensure that it was between 0 and 360.
    !    Then mod(A,360) would do, if A was positive, but if A
    !    was negative, your result would be between -360 and 0.
    !
    !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
    !
    !    An I4 is an integer ( kind = 4 ) value.
    !
    !  Example:
    !
    !        I     J     MOD I4_MODP    Factorization
    !
    !      107    50       7       7    107 =  2 *  50 + 7
    !      107   -50       7       7    107 = -2 * -50 + 7
    !     -107    50      -7      43   -107 = -3 *  50 + 43
    !     -107   -50      -7      43   -107 =  3 * -50 + 43
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the number to be divided.
    !
    !    Input, integer ( kind = 4 ) J, the number that divides I.
    !
    !    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
    !    divided by J.
    !
      implicit none

      integer ( kind = 4 ) i
      integer ( kind = 4 ) i4_modp
      integer ( kind = 4 ) j
      integer ( kind = 4 ) value

      if ( j == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
      end if

      value = mod ( i, j )

      if ( value < 0 ) then
        value = value + abs ( j )
      end if

      i4_modp = value

      return
    end


    function i4_wrap ( ival, ilo, ihi )
    !! I4_WRAP forces an I4 to lie between given limits by wrapping.
    !
    !  Discussion:
    !
    !    An I4 is an integer ( kind = 4 ) value.
    !
    !    There appears to be a bug in the GFORTRAN compiler which can lead to
    !    erroneous results when the first argument of I4_WRAP is an expression.
    !    In particular:
    !
    !    do i = 1, 3
    !      if ( test ) then
    !        i4 = i4_wrap ( i + 1, 1, 3 )
    !      end if
    !    end do
    !
    !    was, when I = 3, returning I4 = 3.  So I had to replace this with
    !
    !    do i = 1, 3
    !      if ( test ) then
    !        i4 = i + 1
    !        i4 = i4_wrap ( i4, 1, 3 )
    !      end if
    !    end do
    !
    !  Example:
    !
    !    ILO = 4, IHI = 8
    !
    !    I  Value
    !
    !    -2     8
    !    -1     4
    !     0     5
    !     1     6
    !     2     7
    !     3     8
    !     4     4
    !     5     5
    !     6     6
    !     7     7
    !     8     8
    !     9     4
    !    10     5
    !    11     6
    !    12     7
    !    13     8
    !    14     4
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    07 September 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) IVAL, a value.
    !
    !    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
    !
    !    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
    !
      implicit none

      integer ( kind = 4 ) i4_wrap
      integer ( kind = 4 ) ihi
      integer ( kind = 4 ) ilo
      integer ( kind = 4 ) ival
      integer ( kind = 4 ) jhi
      integer ( kind = 4 ) jlo
      integer ( kind = 4 ) value
      integer ( kind = 4 ) wide

      jlo = min ( ilo, ihi )
      jhi = max ( ilo, ihi )

      wide = jhi - jlo + 1

      if ( wide == 1 ) then
        value = jlo
      else
        value = jlo + i4_modp ( ival - jlo, wide )
      end if

      i4_wrap = value

      return
    end


    subroutine CalcCoeff(NumDim, N_Idx, Level_Max, N_Forward, Idx, coef)
    !
    !   Calculate the combinatorial coefficient
    !
    !   Args:
    !       NumDim      : Total number of dimensions
    !       N_Idx       : Total numner of incices
    !       Level_Max   : Maximum allowed level
    !       N_Forward   : Forward neighbours
    !       Idx         : The actual indexes
    !   Return:
    !       
    !       coef        : The coefficient for each index
    !
        implicit none

        integer     :: NumDim, N_Idx 
        integer     :: N_Forward(NumDim,N_Idx), Idx(NumDim,N_Idx), coef(N_Idx)
        
        integer     :: nbb
        integer     :: pos, d1, D, Nalpha
        integer     :: alpha(NumDim)
        integer     :: Level_Max
        
        logical     :: neighbour, more
        
        coef(1:N_Idx) = 0
        
        do pos=1, N_Idx

            alpha = 0
            more=.true.
            do while (more)
                call Next_alpha( NumDim, Idx(:,pos), Level_Max, alpha, more )
                
                ! Walk through the "Neighbourhood"
                neighbour = .true.
                nbb=pos
                do d1=1,NumDim
                    if(alpha(d1) .ne. 0) nbb = N_Forward(d1,nbb)
                    if ( nbb < 1 )  then
                        neighbour = .false.
                        exit
                    end if
                end do
                
                if (.not. neighbour) cycle
                
                ! Add (-1)^D to the coefficient
                D=sum( alpha )
                coef(pos) = coef(pos) + 1 - 2 * mod(D,2)
                
            
            end do

        end do


    end subroutine CalcCoeff


    subroutine Next_alpha( NumDim, Idx_pos, Level_Max, alpha, more )
    !
    !   Calculate the next binary vector alpha = (0,0,1,1, ....)
    !
    !
    !   Args:
    !       NumDim
    !       Active_Max
    !   Return:
    !       N_Forward   : Forward neighbours
    !       N_Backward  : Backward neighbours
    !       N_Idx       : Total number of indexes
    !       Idx         : The actual indexes
    !       Indicator   : ==1:Active, ==0, Old
      integer  :: NumDim
      integer  :: N_Idx
      integer  :: pos
      integer  :: Idx_pos(NumDim)
      logical  :: more
      integer  :: alpha(NumDim)
      integer  :: Level_Max
      
      integer  alpha_sum
      integer  i,j
      
      

        do
            i = 0
            do while ( i < NumDim )
                i = i + 1
                if ( alpha(i) == 1 ) then
                    alpha(i) = 0
                else 
                    alpha(i) = 1
                    do
                        if ( sum (Idx_pos + alpha ) <= Level_Max ) then          
                            exit
                        end if
                        
    ! Binary addition
                        alpha(i) = 0
                        do while ( i < NumDim )
                            i = i + 1
                            if ( alpha(i) == 1 ) then
                                alpha(i) = 0
                            else
                                alpha(i) = 1
                                exit
                            end if
                        end do
                    end do
                    exit
                end if
            end do

            
            alpha_sum = sum ( alpha(1:NumDim) )
            if ( alpha_sum == 0 ) then
            more = .false.
                exit
            end if
            
            return

      end do

      return

    end subroutine Next_alpha

end module GerstnerGriebel
