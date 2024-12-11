program main
    
use blas95
use lapack95
use flib_dom

use XMLModelBuilder

use Cable_mod
use Domain_mod
    
implicit none

character(80) :: inpFilename
integer :: cmdLineStatus

call getarg(1, inpFilename, cmdLineStatus)

call buildModel(inpFileName)

call destroy_domain

end program main