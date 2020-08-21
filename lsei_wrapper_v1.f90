subroutine lsei_wrapper(w, wdim, mdw, me, ma, mg, n, prgopt, x, rnorme, rnorml, mode, ws, wsdim, ip) bind(C)

    use iso_c_binding
    implicit none

    integer (c_long), value, intent(in)                 :: mdw, me, ma, mg, n, wdim, wsdim
    real(c_double), dimension(mdw,n+1), intent(inout)   :: w
    real(c_double), intent(out)                         :: rnorme, rnorml
    real(c_double), dimension(30), intent(in)           :: prgopt
    real(c_double), dimension(n+10), intent(out)        :: x
    integer (c_long), dimension(3), intent(inout)       :: ip
    integer (c_long), intent(out)                       :: mode
    real(c_double), dimension(wsdim), intent(out)       :: ws 
    integer :: mdw1, me1, ma1, mg1, n1, mode1, wsdim1
    integer, dimension(3) :: ip1

    interface
          SUBROUTINE LSEI(W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME, RNORML, MODE, WS, IP)
            INTEGER MDW, ME, MA, MG, N, MODE, IP(3)
            DOUBLE PRECISION W(:,:), PRGOPT(30), X(:), WS(:), RNORME, RNORML
          end SUBROUTINE LSEI
    end interface

    ! convert from long integer to short 
    mdw1 = mdw
    me1 = me
    ma1 = ma
    mg1 = mg
    n1 = n
    mode1 = mode
    wsdim1 = wsdim
    ip1(1) = ip(1)
    ip1(2) = ip(2)
    ip1(3) = ip(3)

    call LSEI(w, mdw1, me1, ma1, mg1, n1, prgopt, x, rnorme, rnorml, mode1, ws, ip1)

end subroutine lsei_wrapper
