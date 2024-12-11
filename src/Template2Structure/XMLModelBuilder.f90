module XMLModelBuilder

use flib_dom

use Domain_mod
use Cable_mod

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

! Rigid offsets
domNodes => getElementsByTagName(doc,"RigidOffset")
print*,'Processing rigid offsets ...'
do i=0,getLength(domNodes)-1
    call processRigidOffset(item(domNodes,i))
enddo

! Cables
domNodes => getElementsByTagName(doc,"Cable")
print*,'Processing cables ...'
do i=0,getLength(domNodes)-1
    call processCable(item(domNodes,i))
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

subroutine processRigidOffset(domNode)

type(fNode), pointer, intent(in) :: domNode
character(10) :: attribString
integer(kind=4) :: ID
real(kind=8) :: x, y, z
type(RigidOffset_t), pointer :: offsetptr

attribString = getAttribute(domNode, 'ID')
ID = jnum(attribString) 
attribString = getAttribute(domNode, 'X')
x = dnum(attribString)
attribString = getAttribute(domNode, 'Y')
y = dnum(attribString)
attribString = getAttribute(domNode, 'Z')
z = dnum(attribString)

allocate(offsetptr)
call make_RigidOffset(offsetptr, ID, x, y, z)
call add_rigidoffset_to_domain(offsetptr)

end subroutine processRigidOffset
    
!-------------------------------------------------

subroutine processCable(domNode)

type(fNode), pointer, intent(in) :: domNode
type(fNodeList), pointer :: domNodes
type(fNode), pointer :: refgeomDOMnode
type(fNode), pointer :: propertiesDOMnode
type(fNode), pointer :: sectionInertiaDOMnode
type(fNode), pointer :: splineparamDOMnode
character(20) :: refgeomFile
character(50) :: attribString
integer(kind=4) :: ID
integer(kind=4), dimension(2) :: nodeIDs, rigidOffsetIDs
logical :: active
real(kind=8), dimension(9) :: RE1, RE2

! Spline parameters
integer :: N ! number of control points
integer :: d ! Degree -- cubic spline
integer :: dbrev ! Degree for twist
integer :: dbar ! Degree for strain projection basis
integer :: Ng ! Number of Gauss points per element

real(kind=8), dimension(101,3) :: refGeom ! Coordinates pf reference geometry
integer :: i

! Material and section properties
real(kind=8) :: rho ! mass per unit length
real(kind=8) :: EI 
real(kind=8) :: EA 
real(kind=8) :: GJ 
real(kind=8) :: betBEND ! damping coefficient for bending
real(kind=8) :: betAX ! damping coefficient for axial
real(kind=8) :: betTOR ! damping coefficient for torsion
real(kind=8) :: alph0 ! mass proportional damping factor
real(kind=8), dimension(3,3) :: II

type(Cable_t), pointer :: cableptr

attribString = getAttribute(domNode, 'ID')
ID = jnum(attribString) 
attribString = getAttribute(domNode, 'Nodes')
call csv2array(2, attribString, nodeIDs)
attribString = getAttribute(domNode, 'RigidOffsets')
call csv2array(2, attribString, rigidoffsetIDs)
attribString = getAttribute(domNode, 'Active')
if (attribString == "true") then
    active = .true.
else
    active = .false.
endif
attribString = getAttribute(domNode, 'Frame1')
call get3x3MatrixFromString(attribString, RE1, .false.)
attribString = getAttribute(domNode, 'Frame2')
call get3x3MatrixFromString(attribString, RE2, .false.)

domNodes => getElementsByTagName(domNode,'ReferenceGeometry')
refgeomDOMnode => item(domNodes,0)
refgeomFile = getAttribute(refgeomDOMnode,'File')
domNodes => getElementsByTagName(domNode,'Properties')
propertiesDOMnode => item(domNodes,0)
attribString = getAttribute(propertiesDOMnode,'rho')
rho = dnum(attribString)
attribString = getAttribute(propertiesDOMnode,'EA')
EA = dnum(attribString)
attribString = getAttribute(propertiesDOMnode,'EI')
EI = dnum(attribString)
attribString = getAttribute(propertiesDOMnode,'GJ')
GJ = dnum(attribString)
attribString = getAttribute(propertiesDOMnode,'SectionMassMomentOfInertia')
call get3x3MatrixFromString(attribString, II, .true.)
domNodes => getElementsByTagName(domNode,'SplineParameters')
splineparamDOMnode => item(domNodes,0)
attribString = getAttribute(splineparamDOMnode,'N')
N = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'d')
d = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'dtwist')
dbrev = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'dstrain')
dbar = jnum(attribString)
attribString = getAttribute(splineparamDOMnode,'NGauss')
Ng = jnum(attribString)

! Read reference geometry from file
open(file=trim(refgeomFile), unit=100, status='old')
do i = 1,101
    read(100,*)refGeom(i,1:3)
enddo
close(unit=100)

allocate(cableptr)
call make_cable(cableptr, ID, nodeIDs, rigidoffsetIDs, active, &
                RE1, RE2, &
                N, d, dbrev, dbar, Ng, refGeom, &
                rho, EI, EA, GJ, betBEND, betAX, betTOR, alph0, II)
call add_element_to_domain(cableptr)

end subroutine processCable
    
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