module XMLModelBuilder

use flib_dom

use Domain_mod

implicit none

interface csv2array
    module procedure csv2IntArray
    module procedure csv2DoubleArray
end interface

private csv2array, csvCountElems

contains
    
!-------------------------------------------------
    
subroutine buildModel(filename)

character*(*), intent(in) :: filename

type(fNode), pointer :: doc
type(fNodeList), pointer :: domNodes

integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 DOMAIN BUILD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*,'Processing XML input file ...'
doc => parsefile(filename)

call make_Domain
call startDomainBuild

! Nodes
domNodes => getElementsByTagName(doc,"Node")
print*,'Processing nodes ...'
do i=0,getLength(domNodes)-1
    call processNode(item(domNodes,i))
enddo

call endDomainBuild

call destroyNode(doc)

end subroutine buildModel

!-------------------------------------------------

subroutine processNode(domNode)

type(fNode), pointer, intent(in) :: domNode
character(10) :: attribString
integer(kind=4) :: ID
real(kind=8) :: x, y, z, RJ(9)
integer :: num_spec_constraints ! number of constraints specified in the input
integer :: constraints(NUM_NODE_DOF)
integer :: i
type(Node_t), pointer :: nodeptr

attribString = getAttribute(domNode, 'ID')
ID = jnum(attribString) ! What if ID_TYPE is integer*8 ?
attribString = getAttribute(domNode, 'X')
x = dnum(attribString)
attribString = getAttribute(domNode, 'Y')
y = dnum(attribString)
attribString = getAttribute(domNode, 'Z')
z = dnum(attribString)
attribString = getAttribute(domNode, 'Frame')
call get3x3MatrixFromString(attribString, RJ, .true.)
attribString = getAttribute(domNode, 'Constraints')
num_spec_constraints = len_trim(attribString)
constraints = 0
do i=1,num_spec_constraints
    if (attribString(i:i) .eq. '1') constraints(i) = 1
enddo
if (num_spec_constraints .lt. NUM_NODE_DOF) then
    constraints(num_spec_constraints+1:NUM_NODE_DOF) = 1
endif

allocate(nodeptr)
call make_Node(nodeptr, ID, x, y, z, RJ, constraints)
call add_node_to_domain(nodeptr)

end subroutine processNode
    
!-------------------------------------------------
! Functions related to parsing CSV strings
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