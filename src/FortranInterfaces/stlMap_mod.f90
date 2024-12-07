module stlMap_mod
    
use, intrinsic :: iso_c_binding

implicit none

interface
    function create_map() result(map_ptr) &
      bind(C, name="CREATE_MAP")
        use, intrinsic :: iso_c_binding
        type(c_ptr) :: map_ptr
    end function create_map
end interface

interface
    subroutine delete_map(map_ptr) &
      bind(C, name="DELETE_MAP")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: map_ptr
    end subroutine delete_map
end interface

interface
    subroutine add_to_map(map_ptr, key, val) &
      bind(C, name="ADD_TO_MAP")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: map_ptr
        integer(kind=4), value, intent(in) :: key
        type(c_ptr), value, intent(in) :: val
    end subroutine add_to_map
end interface

interface
    subroutine remove_from_map(map_ptr, key) &
      bind(C, name="ADD_TO_MAP")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: map_ptr
        integer(kind=4), value, intent(in) :: key
    end subroutine remove_from_map
end interface

interface
    function get_val_for_key(map_ptr, key) result(valptr) &
      bind(C, name="GET_VAL_FOR_KEY")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: map_ptr
        integer(kind=4), value, intent(in) :: key
        type(c_ptr) :: valptr
    end function get_val_for_key
end interface

interface
    function get_begin_iterator(map_ptr) result(iter_ptr) &
      bind(C, name="GET_BEGIN_ITERATOR")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: map_ptr
        type(c_ptr) :: iter_ptr
    end function get_begin_iterator
end interface

interface
    subroutine iterate_next(iter_ptr) &
      bind(C, name="ITERATE_NEXT")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: iter_ptr
    end subroutine iterate_next
end interface

interface
    function get_value_for_iterator(iter_ptr) result(valptr) &
      bind(C, name="GET_VALUE_FOR_ITERATOR")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: iter_ptr
        type(c_ptr) :: valptr
    end function get_value_for_iterator
end interface

interface
    subroutine delete_iterator(iter_ptr) &
      bind(C, name="DELETE_ITERATOR")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: iter_ptr
    end subroutine delete_iterator
end interface

end module stlMap_mod
