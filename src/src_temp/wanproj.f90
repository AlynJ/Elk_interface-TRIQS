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
integer ik,ist,istmin,istmax,ispn,ias,nthd,iorb,is,l
integer lm,rlm,lmmax,isym,jsym,jspn,lspl,lm1,lm2,ia,ja,ka
complex(8), allocatable :: symlm(:,:,:,:), a(:,:)
character(50) exfmt,fname,fext
integer ll,lmax, m, rlmmax,lsplsite
complex(8), allocatable :: z(:),zz(:),z1(:,:)
logical sym
logical, allocatable :: done(:)

!find symmetry operation for local to global transformation of projectors
allocate(symlm(ld,ld,norb,nproj))
do iorb=1,norb
! general projector info
  allocate(done(natmtot))
  done(:)=.false.
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!determine the symmetry operation required to transform
!the equivalent atom ia to ja. 
    lm1=(l+1)**2
!set up identity matrix
    allocate(a(lm1,lm1))
    a(:,:)=0.d0
    do lm=1,lm1
      a(lm,lm)=1.d0
    end do
    do ja=1,natoms(is)
! cycle if not equivalent or symmetry has already been found
      if((done(ja)).or.(.not.eqatoms(ia,ja,is))) cycle
! find first symmetry matrix which transforms ia to ja.
      allocate(z1(lm1,lm1))
      do isym=1,nsymcrys
! index to spatial rotation lattice symmetry
        lspl=lsplsymc(isym)
! check that the crystal symmetry is the same as the site
! symmetry
        ka=ieqatom(ia,is,isym)
! want to find a proper symmetric matrix (i.e. does not 
! include inversion symmetry) which transfroms atom ia
! to ja. This is for local <-> global coordinate transformation. 
        if((ka.eq.ja).and.(symlatd(lspl).gt.0.d0)) then
! calculate the symmerty matrix in lm basis
          call rotzflm(symlatc(:,:,lspl),l,l,lm1,lm1,lm1,a,z1)
! Keep desired array subset 
          symlm(1:lmmax,1:lmmax,iorb,ja)=z1(1:lmmax,1:lmmax)
! found symmetry, skip for this index
          done(ja)=.true.
          !write(*,*) symlm
! Exiting the loop
          exit
        endif
      enddo
      deallocate(z1)
    enddo
    deallocate(a)
  enddo
  deallocate(done)
enddo
! begin parallel loop over k-points
! synchronise MPI processes
! TODO: Need to test the mpi parallelisation
!call mpi_barrier(mpicom,ierror)
! loop over reduced k-point set
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  !if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
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
!call mpi_barrier(mpicom,ierror)
deallocate(symlm)
return
end subroutine
!EOC
