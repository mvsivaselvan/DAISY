subroutine SPECIFIC_DESTROY_MAP(map$, numComponents)

type(c_ptr), intent(in) :: map$
integer(kind=4), intent(in) :: numComponents

type(c_ptr) :: iterator$
type(SPECIFIC_COMPONENT_TYPE), pointer :: obj

integer(kind=4) :: i

if (numComponents .gt. 0) then
    iterator$ = get_begin_iterator(map$)
    do i=1,numComponents
        call c_f_pointer(get_value_for_iterator(iterator$), obj)
        call SPECIFIC_DESTROY_ROUTINE(obj)
        deallocate(obj)
        call iterate_next(iterator$)
    enddo
    call delete_iterator(iterator$)
endif
call delete_map(map$)

end subroutine SPECIFIC_DESTROY_MAP