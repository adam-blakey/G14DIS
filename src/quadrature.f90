module quadrature
    use parameters
    implicit none

    !interface functionArguments
    !    subroutine singleArguementDP(real(dp) x)
    !        use parameters
    !    end subroutine
    !end interface

    contains
        !function GaussLegendreQuadrature(f, n)
        function GaussLegendreQuadrature(f_ptr, n)
            integer, intent(in) :: n
            real(dp) :: GaussLegendreQuadrature
            !procedure(singleArguementDP) :: f_ptr#

            abstract interface
                function func(x)
                    import
                    real(dp) :: func
                    real(dp), intent(in) :: x
                end function
            end interface

            procedure(func), pointer, intent(in) :: f_ptr

            real(dp) :: weight
            integer :: i
            real(dp) :: x
            real(dp) :: h

            h = real(2)/n
            GaussLegendreQuadrature = 0

            do i = 1, n
                x = LegendrePolynomialRoot(n, i)
                weight = real(2)/((1-x**2)*(LegendrePolynomialDerivative(x, n))**2)
                GaussLegendreQuadrature = GaussLegendreQuadrature + weight * f_ptr(x)
            end do
        end function GaussLegendreQuadrature

        recursive function LegendrePolynomial(x, n) result(LPoly)
            integer, intent(in) :: n
            real(dp), intent(in) :: x

            real(dp) :: LPoly

            if (n==0) then
                LPoly = 1
            else if (n==1) then
                LPoly = x
            else
                LPoly = real(2*n - 1)/n * x * LegendrePolynomial(x, n-1) - real(n - 1)/n * LegendrePolynomial(x, n-2) 
            end if
        end function LegendrePolynomial

        recursive function LegendrePolynomialDerivative(x, n) result(LPoly)
            integer, intent(in) :: n
            real(dp), intent(in) :: x

            real(dp) :: LPoly

            if (n==0) then
                LPoly = 0
            else if (n==1) then
                LPoly = 1
            else
                LPoly = real(2*n - 1)/(n - 1) * x * LegendrePolynomialDerivative(x, n-1) &
                        - real(n)/(n - 1) * LegendrePolynomialDerivative(x, n-2)
            end if
        end function LegendrePolynomialDerivative

        function LegendrePolynomialRoot(n, i)
            integer, intent(in) :: n
            integer, intent(in) :: i
            real(dp) :: x

            real(dp) :: LegendrePolynomialRoot

            x = -cos((real(2*i - 1)/(2*n))*PI)

            do while (abs(LegendrePolynomial(x, n)) >= 1e-5)
                x = x - LegendrePolynomial(x, n)/LegendrePolynomialDerivative(x, n)
            end do

            LegendrePolynomialRoot = x
        end function LegendrePolynomialRoot

end module quadrature