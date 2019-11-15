module testFunctions
    implicit none

    contains
        function Adam(x)
            use parameters

            real(dp), intent(in) :: x
            real(dp) :: Adam

            Adam = x**2
        end function
end module

program test
    use quadrature
    use linearSystems
    use parameters
    use testFunctions
    use finiteDifference
    implicit none

    abstract interface
        function func(z)
            import
            real(dp) :: func
            real(dp), intent(in) :: z
        end function
    end interface





    procedure(func), pointer :: f_ptr => null()

    f_ptr => Adam

    
    
end program test