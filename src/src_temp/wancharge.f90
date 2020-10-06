!Written by A. D. N. James

subroutine wancharge(tspndg,tlmdg,gf,nproj,projst,ld,lmmax2,nst,idx,subulm,wantrue,w,dmat)
!Calculate and output the charge of the Wannier orbital iorb
use modmain
use moddftu

implicit none
! input
logical, intent(in) :: tspndg,tlmdg,gf
integer, intent(in) :: nproj
integer, intent(in) :: projst(nkpt)
integer, intent(in) :: ld
integer, intent(in) :: lmmax2
integer, intent(in) :: nst
integer, intent(in) :: idx(nstsv,nkpt)
complex(8), intent(in) :: subulm(ld,ld,norb,nproj)
complex(8), intent(in) :: wantrue(ld,nst,nspinor,norb,nproj,nkpt)
complex(8), intent(in) :: w
complex(8), intent(inout) :: dmat(lmmax2,nspinor,lmmax2,nspinor,natmtot)
! local variables
integer ist,jst,ispn,jspn,l,lm1,lm2,ik,i,j,iorb
integer ind1,ind2,lmmax,ia,ja,is,ias,jas,mst
complex(8), allocatable :: z1(:,:), z2(:,:)
complex(8), allocatable :: z3(:,:)
complex(8) chg,z
complex(8) su2(2,2),b(2,2),c(2,2)
real(8) f
logical tsqaz
real(8) th,v1(3),v2(3),v3(3),t1
character(50) exfmt

!loop over each k-point
do iorb=1,norb
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  i=l**2+1
  j=(l+1)**2
  do ia=1,natoms(ia)
    ias=idxas(ia,is)
    do ik=1,nkpt
      mst=projst(ik)
      allocate(z1(lmmax,mst),z2(lmmax,mst),z3(lmmax,lmmax))
      do ispn=1,nspinor
        do jspn=1,nspinor
          if (tspndg.and.(ispn.ne.jspn)) cycle
!set KS state integers into temp complex arrays
          z1(:,:)=wantrue(1:lmmax,1:mst,ispn,iorb,ia,ik)
          z2(:,:)=wantrue(1:lmmax,1:mst,jspn,iorb,ia,ik)
          z3(:,:)=0.d0
!calculate the local Green's function for frequency w
          if(gf) then
            call gwan(lmmax,mst,ik,idx,w,z1,z2,z3)
          else
!calculate wannier charge for ik
            call wanchg(tlmdg,lmmax,mst,ik,idx,z1,z2,z3)
          endif
          dmat(i:j,ispn,i:j,jspn,ias)=dmat(i:j,ispn,i:j,jspn,ias)+z3(:,:)*wkpt(ik)
        end do
      end do
      deallocate(z1,z2,z3)
    end do
  end do
enddo
! symmetrise the wannier density matrices
call symdmat(lmaxdm,lmmax2,dmat)
! output charge in irreducible form
if (lmirep) then
  do iorb=1,norb
    is=orb(iorb,1)
    l=orb(iorb,2)
    lmmax=2*l+1
    i=l**2+1
    j=(l+1)**2
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      allocate(z1(lmmax,lmmax),z2(lmmax,lmmax),z3(lmmax,lmmax))
      do ispn=1,nspinor
        do jspn=1,nspinor
          z1(:,:)=0.d0
          z2(:,:)=subulm(1:lmmax,1:lmmax,iorb,ia)
          z3(:,:)=dmat(i:j,ispn,i:j,jspn,ias)
          call zgemm('N','N',lmmax,lmmax,lmmax,zone,z2,lmmax, &
                 z3,lmmax,zzero,z1,lmmax)
          call zgemm('N','C',lmmax,lmmax,lmmax,zone,z1,lmmax,z2, &
                 lmmax,zzero,z3,lmmax)
          dmat(i:j,ispn,i:j,jspn,ias)=z3(:,:)
        end do
      end do
      deallocate(z1,z2,z3)
    enddo
  enddo
endif
return

contains

subroutine wanchg(tlmdg,lmmax,mst,ik,idx,z1,z2,z3)

use modmain
use moddftu

implicit none
! input
logical, intent(in) :: tlmdg
integer, intent(in) :: lmmax
integer, intent(in) :: mst
integer, intent(in) :: ik
integer, intent(in) :: idx(nstsv,nkpt)
complex(8), intent(in) :: z1(lmmax,mst),z2(lmmax,mst)
complex(8), intent(inout) :: z3(lmmax,mst)
!local
integer lm1,lm2,ist,jst
real(8) f

do lm1=1,lmmax
  do lm2=1,lmmax
    if (tlmdg.and.(lm1.ne.lm2)) cycle
    do ist=1,mst
      jst=idx(ist,ik)
      f=occsv(jst,ik)
!Calculate the Wannier charge including the occupations
      z3(lm1,lm2)=z3(lm1,lm2)+f*z1(lm1,ist)*conjg(z2(lm2,ist))
    end do
  end do
end do
return
end subroutine

subroutine gwan(lmmax,mst,ik,idx,w,z1,z2,z3)

use modmain
use moddftu

implicit none
! input
integer, intent(in) :: lmmax
integer, intent(in) :: mst
integer, intent(in) :: ik
integer, intent(in) :: idx(nstsv,nkpt)
complex(8), intent(in) :: w
complex(8), intent(in) :: z1(lmmax,mst),z2(lmmax,mst)
complex(8), intent(inout) :: z3(lmmax,lmmax)
!local
integer lm1,lm2,ist,jst
real(8) fi,eig
complex(8), allocatable :: gs(:,:),z(:,:)

allocate(gs(mst,mst),z(lmmax,mst))
gs(:,:)=0.d0
!set KS state integers into temp complex arrays
do ist=1,mst
  jst=idx(ist,ik)
  eig=evalsv(jst,ik)-efermi
!KS greens function for ist band
  gs(ist,ist)=occmax/cmplx(dble(w)-eig,aimag(w),8)
end do
!Project the Green's function
call zgemm('N','N',lmmax,mst,mst,zone, &
       z1,lmmax,gs,mst,zzero,z,lmmax)
call zgemm('N','C',lmmax,lmmax,mst,zone, &
       z,lmmax,z2,lmmax,zzero,z3,lmmax)
deallocate(gs,z)

return
end subroutine


end subroutine

