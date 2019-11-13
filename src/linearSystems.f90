module linearSystems
    use parameters
    implicit none

    contains
        function Thomas(a, b, c, d, n)
            integer,                  intent(in) :: n
            real(dp), dimension(n-1), intent(in) :: a
            real(dp), dimension(n),   intent(in) :: b
            real(dp), dimension(n-1), intent(in) :: c
            real(dp), dimension(n),   intent(in) :: d

            real(dp), dimension(n) :: Thomas

            real(dp), dimension(n-1) :: c_
            real(dp), dimension(n)   :: d_

            integer :: i

            c_(1) = c(1)/b(1)
            do i = 2, n-1
                c_(i) = c(i)/(b(i) - c_(i-1)*a(i-1))
            end do

            d_(1) = d(1)/b(1)
            do i = 2, n
                d_(i) = (d(i) - d_(i-1)*a(i-1))/(b(i) - c_(i-1)*a(i-1))
            end do

            Thomas(n) = d_(n)
            do i = n-1, 1, -1
                Thomas(i) = d_(i) - c_(i)*Thomas(i+1)
            end do

        end function Thomas

end module linearSystems