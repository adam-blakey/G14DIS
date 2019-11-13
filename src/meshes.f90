module MMeshes
    use MMesh
    use parameters
    implicit none

    type TMeshes
        integer                 :: count
        TMesh, dimension(count) :: meshList
    end type TMeshes
end module MMeshes