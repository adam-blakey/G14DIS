module parameters
    implicit none

    integer, parameter :: dp = selected_real_kind(15)
    real(dp), parameter :: PI = 3.1415926536
end module parameters

program test
    use parameters
    implicit none

    print *, PI
end program test