      ! --------------------------------------------------
      SUBROUTINE holderkern (qval, Kmat, Umat, Dvec, Ksum, Bmat, &
         & nobs, y, nlam, ulam, eps, &
          & maxit, anlam, npass, jerr, alpmat)
      ! --------------------------------------------------
            IMPLICIT NONE
            ! - - - arg types - - -
            INTEGER :: nobs
            INTEGER :: nlam
            INTEGER :: anlam
            INTEGER :: jerr
            INTEGER :: maxit
            INTEGER :: npass (nlam)
            DOUBLE PRECISION :: eps
            DOUBLE PRECISION :: qval
            DOUBLE PRECISION :: Kmat (nobs, nobs)
            DOUBLE PRECISION :: Umat (nobs, nobs)
            DOUBLE PRECISION :: Dvec (nobs)
            DOUBLE PRECISION :: Ksum (nobs)
            DOUBLE PRECISION :: Bmat (nobs, nobs)
            DOUBLE PRECISION :: y (nobs)
            DOUBLE PRECISION :: ulam (nlam)
            DOUBLE PRECISION :: alpmat (nobs+1, nlam)
            ! - - - local declarations - - -
            INTEGER :: j
            INTEGER :: l
            DOUBLE PRECISION :: K0 (nobs+1, nobs+1)
            DOUBLE PRECISION :: Ki (nobs+1, nobs)
            DOUBLE PRECISION :: yi (nobs+1)
            DOUBLE PRECISION :: mbd
            DOUBLE PRECISION :: minv
            DOUBLE PRECISION :: decib
            DOUBLE PRECISION :: fdr
            DOUBLE PRECISION :: Ginv
            DOUBLE PRECISION :: r (nobs)
            DOUBLE PRECISION :: phi (nobs)
            DOUBLE PRECISION :: alpvec (nobs+1)
            DOUBLE PRECISION :: oalpvec (nobs+1)
            DOUBLE PRECISION :: dif (nobs+1)
            DOUBLE PRECISION :: Ddiag (nobs, nobs)
            DOUBLE PRECISION :: Ainv (nobs, nobs)
            DOUBLE PRECISION :: BAmat (nobs, nobs)
            DOUBLE PRECISION :: Qinv (nobs, nobs)
            DOUBLE PRECISION :: QKsum (nobs)
            DOUBLE PRECISION :: KUinv (nobs+1, nobs+1)
          
            ! - - - begin - - -
            K0 = 0.0D0
            K0(2:(nobs+1),2:(nobs+1)) = Kmat
            Ki = 1.0D0
            Ki(2:(nobs+1),:) = Kmat
            yi = 1.0D0
            yi(2:(nobs+1)) = y
            mbd = (qval+1.0D0)*(qval+1.0D0)/qval
            minv = 1 / mbd
            decib = qval / (qval + 1.0D0)
            fdr = - decib ** (qval + 1.0D0)
            npass = 0
            r = 0.0D0
            alpmat = 0.0D0
            alpvec = 0.0D0
            lambda_loop: DO l = 1,nlam
               dif = 0.0D0
            ! - - - computing Ku inverse - - - 
                Ddiag = 0.0D0
                DO j = 1, nobs
                   Ddiag(j, j) = 1/(Dvec(j)*Dvec(j) &
                   & + 2.0D0*nobs*ulam(l)*Dvec(j)*minv)
                ENDDO
               Ainv = MATMUL(MATMUL(Umat,Ddiag), TRANSPOSE(Umat))
               BAmat = MATMUL(Bmat, Ainv)
               Ginv = 1.0D0
               DO j = 1, nobs
                   Ginv = Ginv + BAmat(j, j)
                ENDDO
               Ginv = 1 / Ginv
               Qinv = Ainv - Ginv * MATMUL(Ainv, BAmat)
               QKsum = MATMUL(Qinv, Ksum) / nobs
               KUinv = 0.0D0
               KUinv(1, 1) = (1.0D0 + DOT_PRODUCT(QKsum, Ksum)) / nobs
               KUinv(1, 2:(nobs+1))  = -QKsum
               KUinv(2:(nobs+1), 1)  = -QKsum
               KUinv(2:(nobs+1), 2:(nobs+1))  = Qinv
               ! update alpha
               update_alpha: DO
                     DO j = 1, nobs
                         IF (r(j) > decib) THEN
                            phi(j) = r(j) ** (- qval - 1) * fdr
                         ELSE
                            phi(j) = -1.0D0
                         END IF
                     ENDDO
                   oalpvec = alpvec
                   alpvec =  oalpvec + (2*nobs*minv) * MATMUL(KUinv, & 
                            & (-ulam(l)*MATMUL(K0, oalpvec) - &
                            & 0.5 * MATMUL(Ki, (y*phi))/nobs))
                   dif = alpvec - oalpvec
                   r = r + y * MATMUL(dif, Ki)
                   IF (sum(dif*dif)/sum(oalpvec*oalpvec) < eps) EXIT
                   npass(l) = npass(l) + 1
                   IF (sum(npass) > maxit) EXIT
               ENDDO update_alpha
               alpmat(:, l) = alpvec
               IF (sum(npass) > maxit) THEN
                  EXIT
                  jerr = -l
               ENDIF
               anlam = l
            ENDDO lambda_loop

      END SUBROUTINE holderkern