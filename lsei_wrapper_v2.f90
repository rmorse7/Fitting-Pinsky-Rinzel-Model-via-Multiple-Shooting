subroutine lsei_wrapper(w, mdw, me, ma, mg, n, prgopt, x, rnorme, rnorml, mode, ws, wsdim, ip, ips) bind(C)

    use iso_c_binding
    implicit none

    integer (c_long), value, intent(in)                 :: mdw, me, ma, mg, n, wsdim, ips
    real(c_double), dimension(mdw,n+1), intent(inout)   :: w
    real(c_double), intent(out)                         :: rnorme, rnorml
    real(c_double), dimension(30), intent(in)           :: prgopt
    real(c_double), dimension(n+10), intent(out)        :: x
    integer (c_long), dimension(ips+10), intent(inout)  :: ip
    integer (c_long), intent(out)                       :: mode
    real(c_double), dimension(wsdim), intent(out)       :: ws 

    interface
          SUBROUTINE LSEI(W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME, RNORML, MODE, WS, IP)
            !INTEGER MDW, ME, MA, MG, N, MODE, IP(3)
            !DOUBLE PRECISION W(:,:), PRGOPT(30), X(:), WS(:), RNORME, RNORML
            use iso_c_binding
            implicit none
            integer (c_long), value, intent(in) :: MDW, ME, MA, MG, N
            integer (c_long), intent(inout) :: mode
            integer (c_long), dimension(:), intent(inout) :: IP
            real (c_double), dimension(:,:), intent(inout) :: W
            real (c_double), intent(out) :: RNORME, RNORML
            real (c_double), dimension(:), intent(out) :: X, WS
            real (c_double), dimension(30), intent(in) :: PRGOPT
          end SUBROUTINE LSEI
    end interface

    call LSEI(w, mdw, me, ma, mg, n, prgopt, x, rnorme, rnorml, mode, ws, ip)

end subroutine lsei_wrapper
