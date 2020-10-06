!Written by A. D. N. James 
subroutine writewanflags(ld,nproj,subulm)
use modmain
implicit none
! arguements
integer, intent(in) :: ld
integer, intent(in) :: nproj
complex(8), intent(inout) :: subulm(ld,ld,norb,nproj)
! local variables
integer so,eqva,iorb,ia,ja,is,l,lmmax,jas,lm
character(50) fname
logical, allocatable :: done(:)

! write PROJ file
fname='PROJ'//trim(filext)
!write to general projector output file
open(805,file=trim(fname),form='FORMATTED')
!output flag for spin-coupled (spin-orbit coupling) calculations
so = 0
if(spcpl) so = 1
write(805,'(I6,I8,I8,I8,I8" : nproj, nkpt, nspinor, spinorb, natmtot")') &
                                            norb,nkpt,nspinor,so,natmtot
do iorb=1,norb
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  write(805,'(I3,I3,I2,I2" : Proj index, Species index, l, lm submatrix size")') &
                                                 iorb,is,l,orb(iorb,3)
!determine the number of equivalent atoms
  allocate(done(natmtot))
  done(:)=.false.
  do ia=1,natoms(is)
    eqva=0
    do ja=1,natoms(is)
      if(eqatoms(ia,ja,is)) eqva=eqva+1
    enddo
    write(805,'(I6," : No. of equivalent atoms")') eqva
    do ja=1,natoms(is)
      if((done(ja)).or.(.not.eqatoms(ia,ja,is))) cycle
! variables for key in PROJ.OUT
      jas=idxas(ja,is)
      write(805,'(I6,I6," : atom, spatom, lm indices")') ja,jas
!Output flags for TRIQS to know how the projectors have been outputted.
!Note that the projectors are output in spherical harmonics
      if(cubic) then
        write(805,'("     1 : Cubic Harmonics")')
      elseif(lmirep) then
        write(805,'("     2 : Irreducible basis")')
        write(805,*)
        do lm=1,lmmax
          write(805,*) dble(subulm(lm,1:lmmax,iorb,ja))
        end do
        write(805,*)
        do lm=1,lmmax
          write(805,*) aimag(subulm(lm,1:lmmax,iorb,ja))
        end do
      else
        write(805,'("     0 : Spherical Harmonics")')
      end if
      write(805,'("")')
      done(ja)=.true.
    enddo
  enddo
  deallocate(done)
enddo
close(805)
return
end subroutine
!EOC
