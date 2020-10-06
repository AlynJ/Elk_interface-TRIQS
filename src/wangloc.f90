!Written by A. D. N. James

subroutine wangloc(tspndg,tlmdg,matsu,nproj,projst,ld,nst,idx,subulm,wantrue,gloc)
!Calculate the local Green's functions and output its Spectral functions

!Note that the performance of this routine is slow due to the symmetrisation of
!the density matrix in wancharge.f90 for each frequency. This can be
!significantly improved when just symmetrising the impurity atom(s) Greem's
!function. This would require an additional subroutine.

use modmain
use moddftu
use modomp
use modmpi

implicit none
! input
logical, intent(in) :: tspndg,tlmdg,matsu
integer, intent(in) :: nproj
integer, intent(in) :: projst(nkpt)
integer, intent(in) :: ld
integer, intent(in) :: nst
integer, intent(in) :: idx(nstsv,nkpt)
complex(8), intent(in) :: subulm(ld,ld,norb,nproj)
complex(8), intent(in) :: wantrue(ld,nst,nspinor,norb,nproj,nkpt)
complex(8), intent(inout) :: gloc(nwplot,ld,nspinor,ld,nspinor,norb,nproj)
! local variables
integer iw,ispn,jspn,iorb,lmi,lmj
integer nthd,lp,is,ia,ias,i,j,l,lmmax,ld2
real(8) dw
complex(8) g,w
complex(8), allocatable :: z1(:,:),z2(:,:),z3(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
real(8), allocatable :: wr(:)

! Calculate the real fequencies 
allocate(wr(nwplot))
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  wr(iw)=dw*dble(iw-1)+wplot(1)
end do

ld2=(lmaxdm+1)**2
allocate(dmat(ld2,nspinor,ld2,nspinor,natmtot))

! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over reduced k-point set
call holdthd(nwplot/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(dmat,w,iorb) &
!$OMP PRIVATE(is,ia,ias,l,lmmax) &
!$OMP PRIVATE(i,j,ispn,jspn) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do iw=1,nwplot
  dmat(:,:,:,:,:)=0.d0
  w=cmplx(wr(iw),swidth,8)
  call wanchg(.false.,.false.,.true.,nproj,projst,ld, &
                 ld2,nst,idx,subulm,wantrue,w,dmat)
  do iorb=1,norb
    is=orb(iorb,1)
    l=orb(iorb,2)
    lmmax=2*l+1
    i=l**2+1
    j=(l+1)**2
    do ia=1,natoms(is)
    !projector properties
      ias=idxas(ia,is)
    !put dmat into local Green's function
      do ispn=1,nspinor
        do jspn=1,nspinor
          gloc(iw,1:lmmax,ispn,1:lmmax,jspn,iorb,ia)=dmat(i:j,ispn,i:j,jspn,ias)
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
deallocate(dmat)
! broadcast local green's function array array to every MPI process
if (np_mpi.gt.1) then
  do iw=1,nwplot
    do lmi=1,ld
      do ispn=1,nspinor
        do lmj=1,ld
          do jspn=1,nspinor
            do iorb=1,norb
              lp=mod(iw-1,np_mpi)
              call mpi_bcast(gloc(iw,lmi,ispn,lmj,jspn,iorb,:),nproj,mpi_double_complex,lp,mpicom,ierror)
            end do
          end do
        end do
      end do
    end do
  end do
end if
!rotate to new basis
if (lmirep) then
  do iw=1,nwplot
    do iorb=1,norb
      is=orb(iorb,1)
      l=orb(iorb,2)
      lmmax=2*l+1
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ispn=1,nspinor
          do jspn=1,nspinor
            allocate(z1(lmmax,lmmax),z2(lmmax,lmmax),z3(lmmax,lmmax))
            z1(:,:)=0.d0
            z2(:,:)=subulm(1:lmmax,1:lmmax,iorb,ia)
            z3(:,:)=gloc(iw,1:lmmax,ispn,1:lmmax,jspn,iorb,ia)
            call zgemm('N','N',lmmax,lmmax,lmmax,zone,z2,lmmax, &
                 z3,lmmax,zzero,z1,lmmax)
            call zgemm('N','C',lmmax,lmmax,lmmax,zone,z1,lmmax,z2, &
                 lmmax,zzero,z3,lmmax)
            gloc(iw,1:lmmax,ispn,1:lmmax,jspn,iorb,ia)=z3(:,:)
            deallocate(z1,z2,z3)
          end do
        end do
      end do
    end do
  end do
endif
deallocate(wr)
return

end subroutine
!EOC
