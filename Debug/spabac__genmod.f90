        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec 02 21:50:16 2019
        MODULE SPABAC__genmod
          INTERFACE 
            SUBROUTINE SPABAC(KV,LOADS,KDIAG)
              REAL(KIND=8), INTENT(IN) :: KV(:)
              REAL(KIND=8), INTENT(INOUT) :: LOADS(0:)
              INTEGER(KIND=4), INTENT(IN) :: KDIAG(:)
            END SUBROUTINE SPABAC
          END INTERFACE 
        END MODULE SPABAC__genmod