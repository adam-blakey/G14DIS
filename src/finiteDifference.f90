!------------------------------------------------------------------------------
! G14DIS, Blakey FEM
!------------------------------------------------------------------------------
!
! MODULE: finiteDifference
! \package athing
!
!> @author
!> <a href="https://adam.blakey.family" target="_blank">Adam Matthew Blakey</a>
!
!> @details
!> Lots of details!
!
! @DESCRIPTION: 
!> This module stores things to do with finite difference schemes.
!
! REVISION HISTORY:
! 3 Nov 2019 - Initial Version
!------------------------------------------------------------------------------

module finiteDifference
    use parameters
    implicit none

    contains
        function poisson1D(ptr_f, N, x0, xN, a, b, c, d)
            abstract interface
                function realDpToRealDp(z)
                    import
                    real(dp) :: func
                    real(dp), intent(in) :: z
                end function
            end interface

            procedure(realDpToRealDp), pointer, intent(in) :: ptr_f
            integer,  intent(in)                           :: N
            real(dp), intent(in)                           :: x0
            real(dp), intent(in)                           :: xN
            real(dp), intent(in), dimension(N-1)           :: a
            real(dp), intent(in), dimension(N)             :: b
            real(dp), intent(in), dimension(N-1)           :: c
            real(dp), intent(in), dimension(N)             :: d

            real(dp), dimension(N+1) :: poisson1D

            real(dp) :: h
            h = (xN - x0)/N

            something(1) = 1

            something(N) = 1


        end function poisson1D

end module finiteDifference