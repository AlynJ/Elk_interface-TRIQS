!Written by A. D. N. James 
subroutine writewanproj(ld,nst,nproj,projst,idx,subulm,sublm,wanprj)
!write the wannier projectors

use modmain
implicit none
!arguments
integer, intent(in) :: ld
integer, intent(in) :: nst
integer, intent(in) :: nproj
integer, intent(in) :: projst(nkpt)
integer, intent(in) :: idx(nstsv,nkpt)
complex(8), intent(in) :: subulm(ld,ld,norb,nproj)
integer, intent(in) :: sublm(ld,norb,nproj)
complex(8), intent(in) :: wanprj(ld,nst,nspinor,norb,nproj,nkpt)
! local variables
integer iorb,ik,l,is,ia,ispn,jspn,lm1,lm2,ias,ja,jas,na
integer mst,ist,istmin,istmax,lm,rlm,lmmax,eqva,so
character(50) fname,extfmt
! allocatable arrays
integer, allocatable :: idxwin(:,:)
logical, allocatable :: done(:)

! write PROJ file
fname='PROJ'//trim(filext)
!write to general projector output file
open(805,file=trim(fname),form='FORMATTED')
!output flag for spin-coupled (spin-orbit coupling) calculations
so = 0
if(spcpl) so = 1
write(805,'(I4,I8,I8,I8,I8" : nproj, nkpt, nspinor, spinorb, natmtot")') &
                                            norb,nkpt,nspinor,so,natmtot
do iorb=1,norb
  is=orb(iorb,1)
  l=orb(iorb,2)
  rlm=orb(iorb,3)
  lmmax=2*l+1
  na=natoms(is)
  write(805,'(I4" : Proj index")') iorb
  write(805,'(I4,I4,I4,I4" : Species index, natoms, l, lm submatrix size")') is,na,l,rlm
!determine the number of equivalent atoms
  allocate(done(natmtot))
  done(:)=.false.
  write(extfmt,'("("I3"(I4)")') rlm
  extfmt=trim(extfmt)//'," : lm indices")'
  do ia=1,natoms(is)
    if(done(ia)) cycle
    eqva=0
    do ja=1,natoms(is)
      if(eqatoms(ia,ja,is)) eqva=eqva+1
    enddo
    write(805,'(I6," : Subset no. of equivalent atoms")') eqva
    do ja=1,natoms(is)
      if(.not.eqatoms(ia,ja,is)) cycle
! variables for key in PROJ.OUT
      jas=idxas(ja,is)
      write(805,'(I4,I4," : atom, spatom")') ja,jas
      write(805,extfmt) sublm(1:rlm,iorb,ia)
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

! write projector files
do iorb=1,norb
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  rlm=orb(iorb,3)
  do ia=1,natoms(is)
    write(fname,'("WANPROJ_L",I2.2,"_S",I2.2,"_A",I4.4)') l,is,ia
    open(800,file=trim(fname)//trim(filext),form='FORMATTED')
    write(800,'(I8,I8,I8," : number of k-points, lmmax, reduced lmmax")') &
                                                             nkpt,lmmax,rlm
    allocate(idxwin(2,nspinor))
    do ik=1,nkpt
      mst=projst(ik)
      idxwin(:,:)=0
! determine the band index window for each spinor
      if((.not.spinorb).and.(spinpol)) then
        do ist=1,mst
          if((idx(ist,ik).gt.nstfv)) then
      !if there are no majority states
            if(ist.eq.1) then
              idxwin(1,2)=1
              idxwin(2,2)=mst
              exit
      !if there are no minority states
            elseif(ist.eq.nst) then
              idxwin(1,1)=1
              idxwin(2,1)=mst
              exit
            else
              idxwin(1,1)=1
              idxwin(2,1)=ist-1
              idxwin(1,2)=ist
              idxwin(2,2)=mst
              exit
            endif
          endif
        enddo
      else
        idxwin(1,:)=1
        idxwin(2,:)=mst
      endif
      write(800,'(I8,3G18.10," : k-point index, k-point (lattice coordinates)")') ik,vkl(:,ik)
      do ispn=1,nspinor
        istmin=idxwin(1,ispn)
        istmax=idxwin(2,ispn)
        if((istmin.eq.0).or.(istmax.eq.0)) then
          write(800,'(I8,I8,I8," : spinor index, minimum and maximum band indices ")') ispn,0,0
          cycle
        endif
    !real part
        write(800,'(I8,I8,I8," : spinor index, minimum and maximum band indices")') &
                                                 ispn,idx(istmin,ik),idx(istmax,ik)
        do lm=1,lmmax
          write(800,*) (dble(wanprj(lm,ist,ispn,iorb,ia,ik)), ist=istmin,istmax)
        end do
        write(800,'("")')
    !imag part
        do lm=1,lmmax
          write(800,*) (aimag(wanprj(lm,ist,ispn,iorb,ia,ik)), ist=istmin,istmax)
        end do
        write(800,'("")')
      end do
    end do
    deallocate(idxwin)
    close(800)
  enddo
enddo
return
end subroutine

