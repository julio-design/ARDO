        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec 02 21:50:15 2019
        MODULE SPABAC_GAUSS__genmod
          INTERFACE 
            SUBROUTINE SPABAC_GAUSS(KV,LOADS,KDIAG)
              REAL(KIND=8), INTENT(IN) :: KV(:)
              REAL(KIND=8), INTENT(INOUT) :: LOADS(0:)
              INTEGER(KIND=4), INTENT(IN) :: KDIAG(:)
            END SUBROUTINE SPABAC_GAUSS
          END INTERFACE 
        END MODULE SPABAC_GAUSS__genmod