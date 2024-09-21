module emxArray
    
use, intrinsic :: iso_c_binding

implicit none
    
type, bind(C) :: emxArray_real_T
    type(c_ptr) :: data_
    type(c_ptr) :: size_
    integer :: allocatedSize
    integer :: numDimensions
    integer :: canFreeData
end type emxArray_real_T

type :: emxArray_1d_wrapper
    real(kind=8), dimension(:), pointer :: data_
    integer, dimension(:), pointer :: size_
    type(emxArray_real_T) :: emx
end type emxArray_1d_wrapper

type :: emxArray_2d_wrapper
    real(kind=8), dimension(:,:), pointer :: data_
    integer, dimension(:), pointer :: size_
    type(emxArray_real_T) :: emx
end type emxArray_2d_wrapper

contains

!----------------------------------------------------

subroutine emxArray_1d_create(this, nrow)

type(emxArray_1d_wrapper), intent(out) :: this
integer, intent(in) :: nrow

allocate(this%data_(nrow))
allocate(this%size_(1))
this%size_(1) = nrow
this%emx%data_ = c_loc(this%data_)
this%emx%size_ = c_loc(this%size_)
this%emx%allocatedSize = nrow
this%emx%numDimensions = 1
this%emx%canFreeData = 0

end subroutine emxArray_1d_create

!----------------------------------------------------

subroutine emxArray_1d_destroy(this)

type(emxArray_1d_wrapper), intent(out) :: this

if(associated(this%data_)) deallocate(this%data_)
if(associated(this%size_)) deallocate(this%size_)

end subroutine emxArray_1d_destroy

!----------------------------------------------------

subroutine emxArray_2d_create(this, nrow, ncol)

type(emxArray_2d_wrapper), intent(out) :: this
integer, intent(in) :: nrow, ncol

allocate(this%data_(nrow,ncol))
allocate(this%size_(2))
this%size_(1) = nrow
this%size_(2) = ncol
this%emx%data_ = c_loc(this%data_)
this%emx%size_ = c_loc(this%size_)
this%emx%allocatedSize = nrow*ncol
this%emx%numDimensions = 2
this%emx%canFreeData = 0

end subroutine emxArray_2d_create

!----------------------------------------------------

subroutine emxArray_2d_destroy(this)

type(emxArray_2d_wrapper), intent(out) :: this

if(associated(this%data_)) deallocate(this%data_)
if(associated(this%size_)) deallocate(this%size_)

end subroutine emxArray_2d_destroy

!----------------------------------------------------
    
end module emxArray