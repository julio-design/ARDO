MODULE rod
        implicit none
        private
        public :: spabac, sparin, fsparv, rod_km, fkdiag, &
                num_to_g, formnf
CONTAINS
        SUBROUTINE spabac(kv,loads,kdiag)
                !
                ! This subroutine performs Cholesky forward and back-substitution
                ! on a symmetric skyline global matrix.
                !
                IMPLICIT NONE
                INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
                REAL(iwp),INTENT(IN)::kv(:)
                REAL(iwp),INTENT(IN OUT)::loads(0:)
                INTEGER,INTENT(IN)::kdiag(:)
                INTEGER::n,i,ki,l,m,j,it,k
                REAL(iwp)::x
                n=UBOUND(kdiag,1)
                loads(1)=loads(1)/kv(1)
                DO i=2,n
                ki=kdiag(i)-i
                l=kdiag(i-1)-ki+1 
                x=loads(i)
                IF(l/=i)THEN
                        m=i-1
                        DO j=l,m 
                        x=x-kv(ki+j)*loads(j)
                        END DO
                END IF
                loads(i)=x/kv(ki+i)
                END DO
                DO it=2,n
                i=n+2-it
                ki=kdiag(i)-i
                x=loads(i)/kv(ki+i)
                loads(i)=x
                l=kdiag(i-1)-ki+1
                IF(l/=i)THEN
                        m=i-1
                        DO k=l,m
                        loads(k)=loads(k)-x*kv(ki+k)
                        END DO
                END IF
                END DO
                loads(1)=loads(1)/kv(1)
                RETURN
        END SUBROUTINE spabac               

        SUBROUTINE sparin(kv,kdiag)
                !
                ! This subroutine performs Cholesky factorisation on a symmetric
                ! skyline global matrix.
                !
                IMPLICIT NONE
                INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
                REAL(iwp),INTENT(IN OUT)::kv(:)
                INTEGER,INTENT(IN)::kdiag(:)
                INTEGER::n,i,ki,l,kj,j,ll,m,k
                REAL(iwp)::x
                n=UBOUND(kdiag,1)  
                kv(1)=SQRT(kv(1))
                DO i=2,n
                ki=kdiag(i)-i
                l=kdiag(i-1)-ki+1
                DO j=l,i
                x=kv(ki+j)
                kj=kdiag(j)-j
                IF(j/=1)THEN
                        ll=kdiag(j-1)-kj+1
                        ll=max(l,ll)
                        IF(ll/=j)THEN
                                m=j-1
                                DO k=ll,m 
                                x=x-kv(ki+k)*kv(kj+k) 
                                END DO
                        END IF
                END IF
                kv(ki+j)=x/kv(kj+j)
                END DO
                kv(ki+i)=SQRT(x)
                END DO
                RETURN
        END SUBROUTINE sparin

        SUBROUTINE fsparv(kv,km,g,kdiag)
                !
                ! This subroutine assembles element matrices into a symmetric skyline
                ! global matrix.
                !
                IMPLICIT NONE
                INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
                INTEGER,INTENT(IN)::g(:),kdiag(:)
                REAL(iwp),INTENT(IN)::km(:,:)
                REAL(iwp),INTENT(OUT)::kv(:) 
                INTEGER::i,idof,k,j,iw,ival
                idof=UBOUND(g,1)
                DO i=1,idof
                k=g(i)
                IF(k/=0)THEN
                        DO j=1,idof
                        IF(g(j)/=0)THEN
                                iw=k-g(j)
                                IF(iw>=0)THEN
                                        ival=kdiag(k)-iw
                                        kv(ival)=kv(ival)+km(i,j) 
                                END IF
                        END IF
                        END DO
                END IF
                END DO
                RETURN
        END SUBROUTINE fsparv

        SUBROUTINE rod_km(km,ea,length)
                !
                ! This subroutine forms the stiffness matrix of a 1-d "rod" element.
                !
                IMPLICIT NONE
                INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
                REAL(iwp),INTENT(IN)::ea,length
                REAL(iwp),INTENT(OUT)::km(:,:)
                REAL(iwp)::one=1.0_iwp
                km(1,1)=one
                km(2,2)=one
                km(1,2)=-one
                km(2,1)=-one
                km=km*ea/length
                RETURN
        END SUBROUTINE rod_km

        SUBROUTINE fkdiag(kdiag,g)
                !
                ! This subroutine computes the skyline profile.
                !
                IMPLICIT NONE
                INTEGER,INTENT(IN)::g(:)
                INTEGER,INTENT(OUT)::kdiag(:)
                INTEGER::idof,i,iwp1,j,im,k
                idof=SIZE(g)
                DO i=1,idof
                iwp1=1
                IF(g(i)/=0)THEN
                        DO j=1,idof
                        IF(g(j)/=0)THEN
                                im=g(i)-g(j)+1
                                IF(im>iwp1)iwp1=im
                        END IF
                        END DO
                        k=g(i)
                        IF(iwp1>kdiag(k))kdiag(k)=iwp1
                END IF
                END DO
                RETURN
        END SUBROUTINE fkdiag

        SUBROUTINE num_to_g(num,nf,g)
                !
                ! This subroutine finds the g vector from num and nf.
                !
                IMPLICIT NONE
                INTEGER,INTENT(IN)::num(:),nf(:,:)  
                INTEGER,INTENT(OUT)::g(:)
                INTEGER::i,k,nod,nodof 
                nod=UBOUND(num,1) 
                nodof=UBOUND(nf,1)
                DO i=1,nod
                k=i*nodof
                g(k-nodof+1:k)=nf(:,num(i))
                END DO
                RETURN
        END SUBROUTINE num_to_g   

        SUBROUTINE formnf(nf)
                ! Used in many subroutines such as p41-p47, p51-p57, p61, p62, p65, p69, p610
                ! p91-p96, p101-p104, p111-p118
                ! This subroutine forms the nf matrix.
                ! nf  : nodal freedom matrix
                IMPLICIT NONE
                INTEGER,INTENT(IN OUT)::nf(:,:)
                INTEGER::i,j,m
                m=0
                DO j=1,UBOUND(nf,2)
                DO i=1,UBOUND(nf,1)
                IF(nf(i,j)/=0)THEN
                        m=m+1
                        nf(i,j)=m
                END IF
                END DO
                END DO
                RETURN
        END SUBROUTINE formnf 
END MODULE rod
