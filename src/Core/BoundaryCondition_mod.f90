module BoundaryCondition_mod

use Node_mod, only : NUM_NODE_DOF, Node_t
    
implicit none
    
type NodeBC_t
    integer(kind=4) :: ID
    integer(kind=4) :: nodeID, toNodeID
    type(Node_t), pointer :: node=>null(), toNode=>null()
    logical, dimension(NUM_NODE_DOF) :: constraints
    real(kind=8), dimension(NUM_NODE_DOF) :: values
    logical :: action ! true = +, false = -
end type NodeBC_t

type BoundaryCondition_t
    integer(kind=4) :: ID
    integer :: numNodeBC = 0
    type(NodeBC_t), dimension(:), pointer :: nodeBCs
    logical :: active
end type BoundaryCondition_t

contains

!--------------------------------------------------------

subroutine make_NodeBC(this, ID, nodeID, constraints, values, toNodeID, action)

type(NodeBC_t), intent(out) :: this
integer(kind=4), intent(in) :: ID, nodeID, toNodeID
logical, dimension(NUM_NODE_DOF), intent(in) :: constraints
real(kind=8), dimension(NUM_NODE_DOF), intent(in) :: values
logical, intent(in) :: action

this%ID = ID
this%nodeID = nodeID
this%toNodeID = toNodeID
this%constraints = constraints
this%values = values
this%action = action

end subroutine make_NodeBC

!--------------------------------------------------------

subroutine destroy_NodeBC(this)

type(NodeBC_t), intent(inout) :: this

if (associated(this%node)) nullify(this%node)
if (associated(this%toNode)) nullify(this%toNode)

end subroutine destroy_NodeBC

!--------------------------------------------------------

subroutine make_BoundaryCondition(this, ID, numNodeBC)

type(BoundaryCondition_t), intent(out) :: this
integer(kind=4), intent(in) :: ID
integer(kind=4), intent(in) :: numNodeBC

this%ID = ID
this%numNodeBC = numNodeBC
allocate( this%nodeBCs(numNodeBC) )
this%active = .false.

end subroutine make_BoundaryCondition

!--------------------------------------------------------

subroutine destroy_BoundaryCondition(this)

type(BoundaryCondition_t), intent(inout) :: this

if (this%numNodeBC .gt. 0) then
    this%numNodeBC = 0
    deallocate(this%nodeBCs)
endif

end subroutine destroy_BoundaryCondition

!--------------------------------------------------------

subroutine activate_BoundaryCondition(this)

use RotationUtils_mod, only : invRodrigues

type(BoundaryCondition_t), intent(inout) :: this
real(kind=8), dimension(NUM_NODE_DOF) :: values

integer(kind=4) :: m, n
type(NodeBC_t), pointer :: nodeBC

this%active = .true.

do m = 1,this%numNodeBC
    nodeBC => this%nodeBCs(m)

    do n = 1,NUM_NODE_DOF
        if (nodeBC%constraints(n)) nodeBC%node%constraints(n) = nodeBC%action
    enddo

    if (.not. nodeBC%action) cycle ! i.e. if removing constraint

    ! The following is if action = .true., i.e. add constraints
    if (associated(nodeBC%toNode)) then ! values are computed from toNode
        values(1:3) = (nodeBC%toNode%x0(1:3) - nodeBC%node%x0(1:3)) &
                    + (nodeBC%toNode%x(1:3) - nodeBC%node%x(1:3))
        values(4:6) = nodeBC%toNode%x(4:6)
    else
        values = nodeBC%values
    endif
    
    do n = 1,NUM_NODE_DOF
        if (nodeBC%constraints(n)) nodeBC%node%x(n) = nodeBC%node%x(n) + values(n)
    enddo
enddo

end subroutine activate_BoundaryCondition

!--------------------------------------------------------

end module BoundaryCondition_mod