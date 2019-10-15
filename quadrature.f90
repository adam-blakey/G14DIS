module quadrature
    use parameters

    recursive function LegendrePolynomial(x, n) result(LPoly)
        integer, parameter:: dp = selected_real_kind(15)

        integer, intent(in):: n
        real(dp), intent(in):: x

        real(dp):: LPoly

        if (n==0) then
            LPoly = 1
        else if (n==1) then
            LPoly = x
        else
            LPoly = real(2*n - 1)/n * x * LegendrePolynomial(x, n-1) - real(n - 1)/n * LegendrePolynomial(x, n-2) 
        end if
    end function LegendrePolynomial

    recursive function LegendrePolynomialDerivative(x, n) result(LPoly)
        integer, parameter:: dp = selected_real_kind(15)
        real(dp), external:: LegendrePolynomial

        integer, intent(in):: n
        real(dp), intent(in):: x

        real(dp):: LPoly

        if (n==0) then
            LPoly = 0
        else if (n==1) then
            LPoly = 1
        else
            LPoly = real(2*n - 1)/(n - 1) * x * LegendrePolynomialDerivative(x, n-1) &
                    - real(n)/(n - 1) * LegendrePolynomialDerivative(x, n-2)
        end if
    end function LegendrePolynomialDerivative

    function LegendrePolynomialRoot(n, i) result(root)
        integer, parameter:: dp = selected_real_kind(15)
        real(dp), parameter:: PI = 3.141592653589793115997963468544185161590576171875

        integer, intent(in):: n
        integer, intent(in):: i
        real(dp):: x

        real(dp), external:: LegendrePolynomial
        real(dp), external:: LegendrePolynomialDerivative

        real(dp):: root

        x = -cos((real(2*i - 1)/(2*n))*PI)

        do while (abs(LegendrePolynomial(x, n)) >= 1e-5)
            x = x - LegendrePolynomial(x, n)/LegendrePolynomialDerivative(x, n)
        end do

        root = x
    end function LegendrePolynomialRoot
end module quadrature