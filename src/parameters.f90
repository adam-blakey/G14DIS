!------------------------------------------------------------------------------
! G14DIS, Blakey FEM
!------------------------------------------------------------------------------
!
! MODULE: Parameters
! \package athing
!
!> @author
!> <a href="https://adam.blakey.family" target="_blank">Adam Matthew Blakey</a>
!
!> @details
!> Lots of details!
!
! @DESCRIPTION: 
!> This module stores various common constant parameters that can be used.
!
! REVISION HISTORY:
! 17 Oct 2019 - Initial Version
! TODO_17_oct_1997 - TODO_YouNeedToDoThings - TODO_Stuff
!------------------------------------------------------------------------------

module parameters
    ! My parameters module.
    implicit none

    !> Double precision constant for 16 bytes of storage.
    integer, parameter :: dp = selected_real_kind(15)

    !> Pi as a constant.
    real(dp), parameter :: PI = 3.1415926536
end module parameters