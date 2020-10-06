subroutine wanchg(tspndg,tlmdg,gf,nproj,projst,ld,ld2,nst,idx,subulm,wantrue,w,dmat)
!Calculate and output the charge of the Wannier orbital iorb
use modmain
use moddftu
!use modomp

implicit none
! input
logical, intent(in) :: tspndg,tlmdg,gf
integer, intent(in) :: nproj
integer, intent(in) :: projst(nkpt)
integer, intent(in) :: ld
integer, intent(in) :: ld2
integer, intent(in) :: nst
integer, intent(in) :: idx(nstsv,nkpt)
complex(8), intent(in) :: subulm(ld,ld,norb,nproj)
complex(8), intent(in) :: wantrue(ld,nst,nspinor,norb,nproj,nkpt)
complex(8), intent(in) :: w
complex(8), intent(inout) ::dmat(ld2,nspinor,ld2,nspinor,natmtot)
! local variables
integer ispn,jspn,l,ik,i,j,iorb
integer lmmax,ia,is,ias,mst
complex(8), allocatable :: z1(:,:), z2(:,:), z3(:,:)
complex(8) chg,z

!loop over each k-point
do ik=1,nkpt
  mst=projst(ik)
  do iorb=1,norb
    is=orb(iorb,1)
    l=orb(iorb,2)
    lmmax=2*l+1
    i=l**2+1
    j=(l+1)**2
    allocate(z1(lmmax,mst),z2(lmmax,mst),z3(lmmax,lmmax))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ispn=1,nspinor
        do jspn=1,nspinor
          if (tspndg.and.(ispn.ne.jspn)) cycle
!set KS state integers into temp complex arrays
          z1(:,:)=wantrue(1:lmmax,1:mst,ispn,iorb,ia,ik)
          z2(:,:)=wantrue(1:lmmax,1:mst,jspn,iorb,ia,ik)
          z3(:,:)=0.d0
!calculate the local Green's function for frequency w
          if(gf) then
            call gwank(lmmax,mst,ik,idx,w,z1,z2,z3)
          else
!calculate wannier charge for ik
            call wanchgk(tlmdg,lmmax,mst,ik,idx,z1,z2,z3)
          endif
          dmat(i:j,ispn,i:j,jspn,ias)=dmat(i:j,ispn,i:j,jspn,ias)+z3(:,:)*wkpt(ik)
        end do
      end do
    end do
    deallocate(z1,z2,z3)
  end do
enddo
! symmetrise the wannier density matrices
call symdmat(lmaxdm,ld2,dmat)
!exit routine for calculating the Green's function
if(gf) return
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
      allocate(z1(lmmax,lmmax),z2(lmmax,lmmax))
      do ispn=1,nspinor
        do jspn=1,nspinor
          z1(:,:)=0.d0
          z2(:,:)=subulm(1:lmmax,1:lmmax,iorb,ia)
          call zgemm('N','N',lmmax,lmmax,lmmax,zone,z2,lmmax, &
                 dmat(i:j,ispn,i:j,jspn,ias),lmmax,zzero,z1,lmmax)
          call zgemm('N','C',lmmax,lmmax,lmmax,zone,z1,lmmax,z2, &
                 lmmax,zzero,dmat(i:j,ispn,i:j,jspn,ias),lmmax)
        end do
      end do
      deallocate(z1,z2)
    end do
  end do
endif
close(810)
return

contains

subroutine wanchgk(tlmdg,lmmax,mst,ik,idx,z1,z2,z3)

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

subroutine gwank(lmmax,mst,ik,idx,w,z1,z2,z3)

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
!EOC
