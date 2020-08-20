subroutine lsei_wrapper(w, mdw, me, ma, mg, n, prgopt, x, rnorme, rnorml, modeLSEI, ws, wsdim, ip, ips) bind(C)

    use iso_c_binding
    implicit none

    integer (c_int32_t), intent(in)                         :: mdw, me, ma, mg, n, wsdim, ips
    real(c_double), dimension(mdw,n+1), intent(inout)       :: w
    real(c_double), intent(out)                             :: rnorme, rnorml
    real(c_double), dimension(30), intent(in)               :: prgopt
    real(c_double), dimension(n+10), intent(out)            :: x
    integer (c_int32_t), dimension(ips+10), intent(inout)   :: ip
    integer (c_int32_t), intent(out)                        :: modeLSEI
    real(c_double), dimension(wsdim), intent(out)           :: ws 

    interface
          SUBROUTINE LSEI(W1, MDW1, ME1, MA1, MG1, N1, PRGOPT1, X1, RNORME1, RNORML1, MODE1, WS1, IP1)
            INTEGER MDW1, ME1, MA1, MG1, N1, MODE1, IP1(:)
            DOUBLE PRECISION W1(:,:), PRGOPT1(30), X1(:), WS1(:), RNORME1, RNORML1
          end SUBROUTINE LSEI
    end interface

    call LSEI(w, mdw, me, ma, mg, n, prgopt, x, rnorme, rnorml, modeLSEI, ws, ip)

end subroutine lsei_wrapper
