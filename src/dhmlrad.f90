
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmlrad
use modmain
use modphonon
implicit none
! local variables
integer is,ias
integer nr,nri,nro,iro,ir
integer l1,l2,l3,m2,lm2
integer npi,i,ilo,jlo,io,jo
real(8) t1,t2
! automatic arrays
real(8) fr1(nrmtmax),fr2(nrmtmax)
! external functions
real(8) splint
external splint
! begin loops over atoms and species
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nro=nr-nri
  iro=nri+1
  npi=npmti(is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
  do l1=0,lmaxapw
    do io=1,apword(l1,is)
      do l3=0,lmaxapw
        do jo=1,apword(l3,is)
          lm2=0
          do l2=0,lmaxi
            do m2=-l2,l2
              lm2=lm2+1
              i=lm2
              do ir=1,nri
                t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*rlmt(ir,2,is)
                fr1(ir)=t1*dble(dvsmt(i,ias))
                fr2(ir)=t1*aimag(dvsmt(i,ias))
                i=i+lmmaxi
              end do
              do ir=iro,nr
                t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*rlmt(ir,2,is)
                fr1(ir)=t1*dble(dvsmt(i,ias))
                fr2(ir)=t1*aimag(dvsmt(i,ias))
                i=i+lmmaxo
              end do
              t1=splint(nr,rlmt(:,1,is),fr1)
              t2=splint(nr,rlmt(:,1,is),fr2)
              dhaa(lm2,jo,l3,io,l1,ias)=cmplx(t1,t2,8)
            end do
          end do
          do l2=lmaxi+1,lmaxo
            do m2=-l2,l2
              lm2=lm2+1
              i=npi+lm2
              do ir=iro,nr
                t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*rlmt(ir,2,is)
                fr1(ir)=t1*dble(dvsmt(i,ias))
                fr2(ir)=t1*aimag(dvsmt(i,ias))
                i=i+lmmaxo
              end do
              t1=splint(nro,rsp(iro,is),fr1(iro))
              t2=splint(nro,rsp(iro,is),fr2(iro))
              dhaa(lm2,jo,l3,io,l1,ias)=cmplx(t1,t2,8)
            end do
          end do
        end do
      end do
    end do
  end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do l3=0,lmaxapw
      do io=1,apword(l3,is)
        lm2=0
        do l2=0,lmaxi
          do m2=-l2,l2
            lm2=lm2+1
            i=lm2
            do ir=1,nri
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*rlmt(ir,2,is)
              fr1(ir)=t1*dble(dvsmt(i,ias))
              fr2(ir)=t1*aimag(dvsmt(i,ias))
              i=i+lmmaxi
            end do
            do ir=iro,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*rlmt(ir,2,is)
              fr1(ir)=t1*dble(dvsmt(i,ias))
              fr2(ir)=t1*aimag(dvsmt(i,ias))
              i=i+lmmaxo
            end do
            t1=splint(nr,rlmt(:,1,is),fr1)
            t2=splint(nr,rlmt(:,1,is),fr2)
            dhloa(lm2,io,l3,ilo,ias)=cmplx(t1,t2,8)
          end do
        end do
        do l2=lmaxi+1,lmaxo
          do m2=-l2,l2
            lm2=lm2+1
            i=npi+lm2
            do ir=iro,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*rlmt(ir,2,is)
              fr1(ir)=t1*dble(dvsmt(i,ias))
              fr2(ir)=t1*aimag(dvsmt(i,ias))
              i=i+lmmaxo
            end do
            t1=splint(nro,rsp(iro,is),fr1(iro))
            t2=splint(nro,rsp(iro,is),fr2(iro))
            dhloa(lm2,io,l3,ilo,ias)=cmplx(t1,t2,8)
          end do
        end do
      end do
    end do
  end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do jlo=1,nlorb(is)
      l3=lorbl(jlo,is)
      lm2=0
      do l2=0,lmaxi
        do m2=-l2,l2
          lm2=lm2+1
          i=lm2
          do ir=1,nri
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*rlmt(ir,2,is)
            fr1(ir)=t1*dble(dvsmt(i,ias))
            fr2(ir)=t1*aimag(dvsmt(i,ias))
            i=i+lmmaxi
          end do
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*rlmt(ir,2,is)
            fr1(ir)=t1*dble(dvsmt(i,ias))
            fr2(ir)=t1*aimag(dvsmt(i,ias))
            i=i+lmmaxo
          end do
          t1=splint(nr,rlmt(:,1,is),fr1)
          t2=splint(nr,rlmt(:,1,is),fr2)
          dhlolo(lm2,jlo,ilo,ias)=cmplx(t1,t2,8)
        end do
      end do
      do l2=lmaxi+1,lmaxo
        do m2=-l2,l2
          lm2=lm2+1
          i=npi+lm2
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*rlmt(ir,2,is)
            fr1(ir)=t1*dble(dvsmt(i,ias))
            fr2(ir)=t1*aimag(dvsmt(i,ias))
            i=i+lmmaxo
          end do
          t1=splint(nro,rsp(iro,is),fr1(iro))
          t2=splint(nro,rsp(iro,is),fr2(iro))
          dhlolo(lm2,jlo,ilo,ias)=cmplx(t1,t2,8)
        end do
      end do
    end do
  end do
! end loops over atoms and species
end do
return
end subroutine

