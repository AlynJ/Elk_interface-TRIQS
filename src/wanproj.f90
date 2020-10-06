!Written by A. D. N. James

subroutine wanproj(ld,nproj,nst,idx,projst,sublm,subulm,wanprj)
!This routine outputs the Wannier projectors for orbital iorb

use modmain
use modomp
use modmpi
implicit none
!arguments
integer, intent(in) :: ld
integer, intent(in) :: nproj
integer, intent(in) :: nst
integer, intent(in) :: idx(nstsv,nkpt)
integer, intent(in) :: projst(nkpt)
integer, intent(in) :: sublm(ld,norb,nproj)
complex(8), intent(in) :: subulm(ld,ld,norb,nproj)
complex(8), intent(inout) :: wanprj(ld,nst,nspinor,norb,nproj,nkpt)
!local
integer i,nea,ik,ist,ispn,ias,jas,lp,nthd,iorb,is,l
integer lm,rlm,lmmax,isym,lspl,lm1,ia,ja,ka,la
complex(8), allocatable :: symlm(:,:,:,:), a(:,:)
complex(8), allocatable :: z(:),zz(:),z1(:,:)
integer, allocatable :: idxiea(:)
logical, allocatable :: done(:)

!find symmetry operation for local to global transformation of projectors
allocate(symlm(ld,ld,norb,nproj))
do iorb=1,norb
! general projector info
  allocate(done(natmtot),idxiea(natmtot))
  done(:)=.false.
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  lm1=(l+1)**2
  allocate(a(lm1,lm1),z1(lm1,lm1))
!set up identity matrix
  a(:,:)=0.d0
  do lm=1,lm1
    a(lm,lm)=1.d0
  end do
!loop over atoms
  do ia=1,natoms(is)
    if(done(ia)) cycle
    ias=idxas(ia,is)
!make an array of equivalent atom indices
    nea=0
    do ja=1,natoms(is)
      if(.not.eqatoms(ia,ja,is)) cycle
      nea=nea+1
      idxiea(nea)=ja
    enddo
    ka=idxiea(1)
!determine the symmetry operation required to transform
!the equivalent atom ja to idxiea(1). 
    do i=1,nea
      ja=idxiea(i)
! find first symmetry matrix which transforms ia to ja.
      do isym=1,nsymcrys
! index to spatial rotation lattice symmetry
        lspl=lsplsymc(isym)
! check that the crystal symmetry is the same as the site
! symmetry
        la=ieqatom(ja,is,isym)
! want to find a proper symmetric matrix (i.e. does not 
! include inversion symmetry) which transfroms atom la
! to ka. This is for local <-> global coordinate transformation. 
        if((ka.eq.la).and.(symlatd(lspl).gt.0.d0)) then
! calculate the symmerty matrix in lm basis
          call rotzflm(symlatc(:,:,lspl),l,l,lm1,lm1,lm1,a,z1)
! Keep desired array subset 
          symlm(1:lmmax,1:lmmax,iorb,ja)=z1(1:lmmax,1:lmmax)
! found symmetry, skip for this index
          done(ja)=.true.
! Exiting the loop
          exit
        endif
      enddo
    enddo
  enddo
  deallocate(done,idxiea,a,z1)
enddo
! begin parallel loop over k-points
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over reduced k-point set
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(wanproj_)
  write(*,'("Info(wanproj): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(wanproj_)
! determine the projectors
  call wanprojk(ik,idx,ld,nproj,nst,projst(ik), &
               subulm,sublm,symlm,wanprj(:,:,:,:,:,ik))
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! broadcast projector array to every MPI process
if (np_mpi.gt.1) then
  do ik=1,nkpt
    do i=1,ld
      do ist=1,nst
        do ispn=1,nspinor
          do iorb=1,norb
            lp=mod(ik-1,np_mpi)
            call mpi_bcast(wanprj(i,ist,ispn,iorb,:,ik),nproj,mpi_double_complex,lp,mpicom,ierror)
          end do
        end do
      end do
    end do
  end do
end if
deallocate(symlm)
return
end subroutine
!EOC
