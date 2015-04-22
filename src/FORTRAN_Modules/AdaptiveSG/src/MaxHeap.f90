module MaxHeap

    implicit none

    contains

    SUBROUTINE heap_adjust(heap, value, root, length)
    !
    !
    !
    !   Args:
    !       root    : The parent node


    real, DIMENSION(:) :: value
    INTEGER, DIMENSION(:) :: heap
    INTEGER :: length, root

    INTEGER :: e, j
    real    :: k

    e=heap(root); k=value(heap(root))

    j=2*root
    DO WHILE (j <= length)
    
       IF(j < length) THEN
          IF(value(heap(j) ) < value(heap(j+1)))j=j+1
       ENDIF
       IF ( k >= value(heap(j)) ) EXIT
       heap(j/2) = heap(j)
       j=j*2
       
    ENDDO
    
    heap(j/2)=e

    END SUBROUTINE heap_adjust


    SUBROUTINE heap_init(heap,value,length, maxvals)

    !
    ! Converts an array into a heap, with the
    ! largest element at the root.
    !

    IMPLICIT NONE

    INTEGER :: length, maxvals
    real :: value(maxvals)
    INTEGER :: heap(length)
    

    INTEGER:: i

    !f2py intent(in)    value, length, maxvals, heap
    !f2py intent(out)  heap

    i=length/2
    DO WHILE(i >= 1)
       CALL heap_adjust(heap,value,i,length)
       i=i-1
    ENDDO

    END SUBROUTINE heap_init
    
    subroutine heap_delete(heap, value, length, rootval)

    !
    ! Removes the root of the heap, then adjusts
    ! the heap. Returns the root contents, i.e. HEAP(1).
    ! Heap cannot be emty, i.e. length > 0 is expected.
    !

    IMPLICIT NONE

    
    INTEGER :: length, maxvals
    real :: value(:)
    INTEGER :: heap(:)
    
    
    
    
    integer :: rootval

    INTEGER :: k, i, j

    !f2py intent(in)    value, length
    !f2py intent(in,out) heap
    !f2py intent(out) rootval

    rootval = heap(1) ! get root contents
    k = heap(length) ! last element must be moved
    length = length - 1 ! reduce heap length by 1

    i = 1; j = 2
    DO WHILE(j.LE.length)
       IF(j.LT.length) THEN
          IF(value(heap(j)) < value(heap(j+1)))j=j+1
       ENDIF  ! j points to the smallest child
       IF(value(k) >= value(heap(j)))EXIT
       heap(i)=heap(j) ! move child up
       i = j; j = 2*j  ! move i and j down
    ENDDO
    heap(i)=k
    

  END subroutine heap_delete
  
  
  SUBROUTINE heap_insert(value, idx, heap, length, maxvals)

    !
    ! Inserts VALUE(idx) in the HEAP.
    ! The length of HEAP is increased by 1.
    ! An empty HEAP, i.e. one of length 0 is allowed.
    !
    ! This has been modified to account for C/python arrays starting
    ! at 0. The value(idx) -> value(idx+1)
    

    IMPLICIT NONE
    
    INTEGER :: idx, length, maxvals
    real :: value(0:maxvals-1)
    INTEGER :: heap(:)
    

    INTEGER :: i, pos
    
    !f2py intent(in)    value, idx, maxvals
    !f2py intent(in,out) heap, length


    length=length+1; i=length
    
    DO WHILE(i.LE.length.AND.i.GE.1)
       IF(i.EQ.1) THEN
          pos=i    ! At root
          EXIT
       ENDIF
       IF(value(idx) <= value(heap(i/2)))THEN
          pos=i
          EXIT
       ENDIF
       heap(i)=heap(i/2) ! Move to parent of i
       i = i/2
    ENDDO

    heap(pos) = idx


    

  END SUBROUTINE heap_insert
  

end module MaxHeap
