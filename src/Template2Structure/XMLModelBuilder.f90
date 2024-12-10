module XMLModelBuilder

implicit none

interface csv2array
    module procedure csv2IntArray
    module procedure csv2DoubleArray
end interface

contains
    
!-------------------------------------------------

subroutine get3x3MatrixFromString(str, mat, sym)

character*(*), intent(in) :: str
real(kind=8), dimension(*), intent(out) :: mat
logical, intent(in) :: sym ! true if symmetric, false if not

integer :: numElems

if (sym) then
    numElems = 6
else
    numElems = 9
endif

call csv2array(numElems, str, mat)

if (sym) then ! fill the upper triangular block
    mat(9) = mat(6)
    mat(6) = mat(5)
    mat(5) = mat(4)
    mat(4) = mat(2)
    mat(7) = mat(3)
    mat(8) = mat(6)
endif

end subroutine get3x3MatrixFromString

!-------------------------------------------------
! Functions related to parsing CSV strings
!-------------------------------------------------

function csvCountElems(csvString) result(numElems)

character*(*), intent(in) :: csvString
integer :: numElems

integer :: lenString, prevIndex

numElems = 0

lenString = len(trim(csvString))
if (lenString .eq. 0) return

prevIndex = lenString+1
do while (prevIndex .gt. 0)
    prevIndex = index(csvString(1:prevIndex-1),',',back=.true.)
    numElems = numElems + 1
enddo

end function csvCountElems

#define SPECIFIC_CSV2ARRAY csv2IntArray
#define SPECIFIC_TYPEKIND integer(kind=4)
#define SPECIFIC_CONVERTER jnum
#include "csv2array_generic.inc"
#undef SPECIFIC_CSV2ARRAY
#undef SPECIFIC_TYPEKIND
#undef SPECIFIC_CONVERTER

#define SPECIFIC_CSV2ARRAY csv2DoubleArray
#define SPECIFIC_TYPEKIND real(kind=8)
#define SPECIFIC_CONVERTER dnum
#include "csv2array_generic.inc"
#undef SPECIFIC_CSV2ARRAY
#undef SPECIFIC_TYPEKIND
#undef SPECIFIC_CONVERTER

!-------------------------------------------------
    
end module XMLModelBuilder