module Reduced_SG

use SandiaRules

implicit none


contains

!#######################################################################################################################################
!
!#######################################################################################################################################

subroutine Max_Next_Points(Idx, N_Idx, dim_num, Coeff, rule, growth, point_total_num)

    integer (kind=4) :: N_Idx, dim_num, point_total_num
    integer (kind=4) :: Idx(dim_num, N_Idx), rule(dim_num), growth(dim_num)
    real (kind=8)    :: Coeff(N_Idx)
    
    integer (kind=4) :: i
    integer (kind=4) :: order_1d(dim_num)
    
    
    !f2py intent(in) Idx, N_Idx, dim_num, Coeff, rule, growth
    !f2py intent(out) point_total_num
    
    point_total_num = 0
    do i=1, N_Idx
        if (Coeff(i) == 0.0) cycle
        call level_growth_to_order ( dim_num, Idx(:,i), rule, growth, order_1d )
        point_total_num = point_total_num + product ( order_1d(1:dim_num) )
    end do
    

end subroutine Max_Next_Points

!#######################################################################################################################################
!
!#######################################################################################################################################

subroutine Calculate_Coefficients(dim_num, Idx, N_Idx, Coeff, q_max)


implicit none

    integer (kind=4)     :: dim_num
    integer (kind=4)     :: N_Idx
    integer (kind=4)     :: Idx(dim_num, N_Idx)
    real (kind=8)        :: Coeff(N_Idx)
    
    integer (kind=4)     :: i
    real (kind=8)        :: q_max
    
    !f2py intent(in) Idx, N_Idx, dim_num, q_max
    !f2py intent(out) Coeff
    
    do i=1, N_Idx
        call sgmga_vcn_coef_DA ( dim_num, Idx, i, N_Idx, q_max, Coeff(i) )
    end do
    
end subroutine Calculate_Coefficients

!#######################################################################################################################################
!
!#######################################################################################################################################

subroutine Weights_and_Points(dim_num, Total_Points, N_Idx, Level_max, Idx, &
Coeff, growth, rule, np, p, Grid_Points)

    integer ( kind=4 )    :: dim_num, Total_Points, Level_max, N_Idx
    real ( kind=8 )       :: Grid_Points(dim_num, Total_Points), Coeff(N_Idx)
    integer ( kind = 4 )  :: Idx(dim_num, N_Idx)
    integer ( kind = 4 )  :: growth(dim_num)
    integer ( kind = 4 )  :: rule(dim_num)
    integer ( kind = 4 ) np(dim_num)
    real ( kind = 8 ) p(*)

  integer ( kind = 4 ) dim
  
  integer     :: level
  logical more_points
  
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_1d(dim_num)
  
  integer ( kind = 4 ) p_index
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_index(dim_num)
  integer ( kind = 4 ) point_total_num2, point_total_num3
  real ( kind = 8 ), allocatable, dimension ( : ) :: points

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_total_index
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_total_order
  
  integer ( kind = 4 ) :: levell(1)
  integer ( kind = 4 ) :: orderl(1)  
  
  integer ( kind = 4 )  :: i

!f2py intent(in) dim_num, Total_Points, N_Idx, Level_max, Coeff, growth, rule, np, p, tol
!f2py intent(out) Grid_Points, Grid_Weights
  

!  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
!  for the TOTAL set of points.
!
  allocate ( sparse_total_order(1:dim_num,1:Total_Points ) )
  allocate ( sparse_total_index(1:dim_num,1:Total_Points ) )

  point_total_num2 = 0
  point_total_num3 = 0
 
  do i=1,N_Idx


    if ( Coeff(i) == 0.0D+00 ) then
      cycle
    end if
!
!  Transform each 1D level to a corresponding 1D order.
!
    call level_growth_to_order ( dim_num, Idx(:,i), rule, growth, order_1d )
    
!
!  The inner loop generates a POINT of the GRID of the LEVEL.
!
    more_points = .false.

    do

      call vec_colex_next3 ( dim_num, order_1d, point_index, more_points )

      if ( .not. more_points ) then
        exit
      end if

      point_total_num2 = point_total_num2 + 1
      sparse_total_order(1:dim_num,point_total_num2) = order_1d(1:dim_num)
      sparse_total_index(1:dim_num,point_total_num2) = point_index(1:dim_num)

    end do
    
!    order_nd = product ( order_1d(1:dim_num) )

!    allocate ( grid_weight(1:order_nd) )

!    call sgmga_product_weight_DA ( dim_num, order_1d, order_nd, rule, &
!      np, p, grid_weight )

!    do order = 1, order_nd
!        point_total_num3 = point_total_num3 + 1
!        Grid_Weights(point_total_num3) = grid_weight(order)
!    end do

!    deallocate ( grid_weight )

  end do
  
  if(point_total_num2 /= Total_Points ) then
    print *, "Total number of points is wrong!"
    stop
  end if 
  
  p_index = 1

  do dim = 1, dim_num

    do level = 0, level_max
	
      levell(1) = level
	  
      call level_growth_to_order ( 1, levell, rule(dim), growth(dim), orderl )
      order=orderl(1)
      allocate ( points(1:order) )

      if ( rule(dim) == 1 ) then
        call clenshaw_curtis_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 2 ) then
        call fejer2_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 3 ) then
        call patterson_lookup_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 4 ) then
         call legendre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 5 ) then
        call hermite_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 6 ) then
        call gen_hermite_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 7 ) then
        call laguerre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 8 ) then
        call gen_laguerre_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 9 ) then
        call jacobi_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 10 ) then
        call hermite_genz_keister_lookup_points_np ( order, np(dim), &
          p(p_index), points )
      else if ( rule(dim) == 11 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 11.'
        stop
      else if ( rule(dim) == 12 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 12.'
        stop
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_SIZE - Fatal error!'
        write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
        stop
      end if

      do point = 1, Total_Points
        if ( sparse_total_order(dim,point) == order ) then
          Grid_Points(dim,point) = &
            points ( sparse_total_index(dim,point) )
        end if
      end do

      deallocate ( points )

    end do

    p_index = p_index + np(dim)

  end do

  deallocate ( sparse_total_index )
  deallocate ( sparse_total_order )

end subroutine Weights_and_Points

!#######################################################################################################################################
!
!#######################################################################################################################################

subroutine Unique_Points(dim_num, Total_Points, seed, tol, Points, u_points, sparse_unique_index)

    integer ( kind = 4 )         :: dim_num, Total_Points, seed
    real (kind=8)                :: Points(dim_num, Total_Points), tol
    integer ( kind = 4 )         :: u_points
    integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
    integer ( kind = 4 )         :: sparse_unique_index(Total_Points)
    
!f2py intent(in) dim_num, Total_Points, Points, tol, seed
!f2py intent(out)   u_points, sparse_unique_index
    
      seed = 123456789
 
  allocate ( undx(1:Total_Points) )

  call point_radial_tol_unique_index ( dim_num, Total_Points, Points, tol, seed, u_points, undx, sparse_unique_index )
  
  deallocate ( undx )
  
end subroutine Unique_Points


!#######################################################################################################################################
!
!#######################################################################################################################################

subroutine Reduce_Points_and_Weights(u_points, dim_num, Total_Points, Points, N_Idx, Idx, &
                                     sparse_unique_index, Coeff, growth, rule, np, p, New_Points, New_Weights)

    integer ( kind = 4 )         :: dim_num, Total_Points, u_points, N_Idx
    real (kind=8)                :: Points(dim_num, Total_Points), Coeff(N_Idx)
    real (kind=8)                :: New_Points(dim_num, u_points), New_Weights(u_points)
    integer ( kind = 4 )         :: sparse_unique_index(Total_Points)
    integer ( kind = 4 )         :: Idx(dim_num, N_Idx), growth(dim_num), rule(dim_num), np(dim_num)
    real (kind=8)                :: p(*)

    integer ( kind = 4 )    :: i, rep
    integer ( kind = 4 )    :: point_total_num
    real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
    real ( kind = 8 )   :: Grid_Weights(Total_Points)
    
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) order_nd

!f2py intent(in) u_points, dim_num, Total_Points, Points, N_Idx, Idx
!f2py intent(in) sparse_unique_index, Coeff, growth, rule, np, p
!f2py intent(out) New_Points, New_Weights

  point_total_num = 0  

  do i=1,N_Idx

    if ( Coeff(i) == 0.0D+00 ) then
      cycle
    end if

    call level_growth_to_order ( dim_num, Idx(:,i), rule, growth, order_1d )
    
    order_nd = product ( order_1d(1:dim_num) )

    allocate ( grid_weight(1:order_nd) )

    call sgmga_product_weight_DA ( dim_num, order_1d, order_nd, rule, &
      np, p, grid_weight )

    do order = 1, order_nd
        point_total_num = point_total_num + 1
        Grid_Weights(point_total_num) = Coeff(i)*grid_weight(order)
    end do

    deallocate ( grid_weight )

  end do


    New_Weights = 0.0d0
    do i=1, Total_Points
        rep = sparse_unique_index(i)
        New_Weights(rep) = New_Weights(rep) + Grid_Weights(i)
        New_Points(:,rep) = Points(:,i)
    end do

end subroutine Reduce_Points_and_Weights

!#######################################################################################################################################
!
!#######################################################################################################################################

!subroutine Concatenate(dim_num, NP1, NP2, Max_Points, Old_Points, New_Points, xdnu2, unique2 )
subroutine Concatenate(dim_num, NP1, NP2, Old_Points, New_Points, xdnu2, unique2 )

    
!    integer ( kind = 4 )         :: dim_num, NP1, NP2, Max_Points
!    real (kind=8)                :: Old_Points(dim_num,Max_Points), New_Points(dim_num, NP2)
    integer ( kind = 4 )         :: dim_num, NP1, NP2
    real (kind=8)                :: Old_Points(dim_num,NP1), New_Points(dim_num, NP2)
    
    integer ( kind = 4 )         :: indx1(NP1), indx2(NP2), indx3(NP1+NP2)
    integer ( kind = 4 )         :: i1, i2, i3
    real (kind=8)                :: r1(NP1), r2(NP2), r3(NP1+NP2)
    integer ( kind = 4 )         :: undx1(NP1), undx2(NP2), undx3(NP1+NP2)
    integer ( kind = 4 )         :: unique_num1, unique_num2, unique_num3
    logical                      :: unique1(NP1), unique2(NP2), unique3(NP1+NP2)
    integer ( kind = 4 )         :: xdnu1(NP1), xdnu2(NP2), xdnu3(NP1+NP2)
    real (kind=8)                :: z(dim_num), a3(dim_num, NP1+NP2)


!f2py intent(in)    dim_num, NP2, Old_Points, NP1, New_Points
!f2py intent(out)   xdnu2, unique2
    
    integer ( kind = 4 ) j
    integer ( kind = 4 ) n3
    integer ( kind = 4 ) seed
    real ( kind = 8 ) tol
    integer ( kind = 4 ) undx_value
    
    seed = 123456789
    tol = sqrt ( r8_epsilon ( ) )
    
    
    call point_radial_tol_unique_index_inc1 ( dim_num, NP1, Old_Points, tol, seed, z, r1, &
    indx1, unique1, unique_num1, undx1, xdnu1 )
    
    do j=1, NP1
        xdnu1(j) = j
    end do
    
    call point_radial_tol_unique_index_inc2 ( dim_num, NP1, Old_Points, NP2, New_Points, tol, z, &
    r1, indx1, unique1, unique_num1, undx1, xdnu1, &
    r2, indx2, unique2, unique_num2, undx2, xdnu2 )

end subroutine Concatenate


!subroutine Calculate_Grid(dim_num, N_Idx, Level_Max, Idx, Coeff, rule, growth, np, p)

!    ! Arguments
!    integer ( kind = 4 ) :: dim_num, N_Idx, Level_Max
!    integer ( kind = 4 ) :: Idx(1:dim_num,1:N_Idx)
!    integer ( kind = 4 ) :: rule(1:dim_num), growth(1:dim_num)
!    integer ( kind = 4 ) :: np(1:dim_num)
!    real ( kind = 8 )    :: p(*), New_Points(:,:), New_Weights(:), Coeff(1:N_Idx)

!    ! Locals
!    integer ( kind = 4 ), parameter    :: seed = 123456789
!    integer ( kind = 4 )               :: point_total_num, unique_point_num
!    real ( kind = 8 ), allocatable     :: Grid_Points(:,:)
!    integer ( kind = 4 ), allocatable  :: sparse_unique_index(:), undx(:)
!    real ( kind = 8 )                  :: tol

!    call Max_Next_Points(Idx, N_Idx, dim_num, Coeff, rule, growth, point_total_num)

!    allocate( Grid_Points(dim_num,point_total_num) )
!    
!    call Weights_and_Points(dim_num, point_total_num, N_Idx, Level_max, Idx, &
!                            Coeff, growth, rule, np, p, Grid_Points)

!    allocate ( undx(1:point_total_num), sparse_unique_index(1:point_total_num) )
!    tol = sqrt ( r8_epsilon ( ) )

!    call point_radial_tol_unique_index ( dim_num, point_total_num, Grid_Points, tol, seed,&
!                                         unique_point_num, undx, sparse_unique_index )
!  
!    deallocate ( undx )
!    
!!    allocate( New_Points(1:dim_num,1:unique_point_num), New_Weights(1:unique_point_num) )
!    
!    call Reduce_Points_and_Weights(unique_point_num, dim_num, point_total_num, Grid_Points, N_Idx, Idx, &
!                                   sparse_unique_index, Coeff, growth, rule, np, p, New_Points, New_Weights)

!    deallocate( Grid_Points, sparse_unique_index )

!end subroutine Calculate_Grid





end module Reduced_SG
