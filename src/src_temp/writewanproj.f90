!Written by A. D. N. James 
subroutine writewanproj(ld,nst,nproj,projst,idx,wanprj)
!write the wannier projectors

use modmain
implicit none
!arguments
integer, intent(in) :: ld
integer, intent(in) :: nst
integer, intent(in) :: nproj
integer, intent(in) :: projst(nkpt)
integer, intent(in) :: idx(nstsv,nkpt)
complex(8), intent(in) :: wanprj(ld,nst,nspinor,norb,nproj,nkpt)
! local variables
integer iorb,ik,l,is,ia,ispn,jspn
integer mst,ist,istmin,istmax,lm,rlm,lmmax
! allocatable arrays
integer, allocatable :: idxwin(:,:)
character(50) fname

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

