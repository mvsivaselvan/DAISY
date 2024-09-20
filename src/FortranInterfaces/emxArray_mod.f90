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
    
end module emxArray