module MMesh
    use parameters
    implicit none

    type TMesh
        integer :: dimension
    end type TMesh

    type, extends(TMesh) :: TInterval
        real(dp) :: start
        real(dp) :: end
    end type TInterval
end module MMesh