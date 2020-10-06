!Written by A. D. N. James 
subroutine writewansym(is,ia,l)
!Write symmetries of each orbital to file.

use modmain
implicit none
! arguments
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: l
! local variables
integer lmi,lmmax,lmmax2,i,j,js,ja,ka
integer isym,lspl,prv_atms
real(8) t0
! allocatable arrays
complex(8), allocatable :: zfmt1(:,:),zfmt2(:,:)
complex(8), allocatable :: symlm(:,:), a(:,:)
integer, allocatable :: ateqv(:,:)
!output string
character(50) fname,exfmt

!atom w.r.t symmetry list
allocate(ateqv(natmtot,nsymcrys))
do isym=1,nsymcrys
! spatial rotation element in lattice point group
  ka=1
  prv_atms=0
  do js=1,nspecies
    do ja=1,natoms(js)
      ateqv(ka,isym)=ieqatom(ja,js,isym)+prv_atms
      ka=ka+1
    end do
    prv_atms=prv_atms+natoms(js)
  end do
end do
!write sym
lmmax=(l+1)**2
lmmax2=2*l+1
lmi=lmmax+1-lmmax2
allocate(a(lmmax,lmmax))
a(:,:)=0.d0
do i=1,lmmax
  a(i,i)=1.d0
end do
! output the l rotation matricies used for projected orbital
allocate(symlm(lmmax,lmmax))
!output symlm for all equivalent atoms 
write(fname,'("PROJSYM_L",I2.2,"_S",I2.2,"_A",I4.4)') l,is,ia
open(820,file=trim(fname)//'.OUT',form='FORMATTED')
write(820,'(I8)') nsymcrys
do isym=1,nsymcrys
  symlm(:,:)=0.d0
! index to spatial rotation lattice symmetry
  lspl=lsplsymc(isym)
! apply the rotation to the muffin-tin function
  call rotzflm(symlatc(:,:,lspl),l,l,lmmax,lmmax,lmmax,a,symlm)
  write(820,'(I8)') isym
  write(820,*) ateqv(:,isym) !atom order wrt symmetry
    !output the lm symmetry matricies
  write(exfmt,'("(",I8,"(G18.10))")') lmmax2
    !real
  do i=1,lmmax2
    write(820,exfmt) (dble(symlm(i,j)), j=1,lmmax2)
  end do
  write(820,'("")')
    !imaginary
  do i=1,lmmax2
    write(820,exfmt) (aimag(symlm(i,j)), j=1,lmmax2)
  end do
  write(820,'("")')
end do
close(820)
deallocate(symlm,a)
return
end subroutine

