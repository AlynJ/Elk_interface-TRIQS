!Written by A. D. N. James

subroutine getdmatdmft

!ADNJ - This routine is the temporary interface from TRIQS to Elk. It reads the DMFT density
! matrix from DMATDMFT.OUT to interface back into the DFT cycle for a self-consistent DFT+DMFT  

use modmain
!use modomp
use modmpi

implicit none
!local variables
integer nkpt_,nspin_, ik, ik_, ispn, ispn_, minst, maxst
integer nst, ist, jst, nspin
real(8) mu, beta
real(8), allocatable :: rd(:), id(:)

!The density matrix is a mixture of the spinor wavefunctions (w.r.t nstsv) for spin-orb,
!therefore this must be accounted for when reading in the density matrix.
if(spinorb) then
  nspin=nspinor-1
else
  nspin=nspinor
endif
! Read the density matrices from DMATDMFT.OUT file
open(32,file='DMATDMFT.OUT',action='READ',form='FORMATTED')
!reads chemical potential (eV)
read(32,*)nkpt_,nspin_,nst,beta,mu
!A check of the read in parameters
if (nkpt.ne.nkpt_) then
  write(*,*)
  write(*,'("Error(getdmatdmft): differing nkpt ")')
  write(*,'(" current    : ",I8)') nkpt
  write(*,'(" DMATDMFT.OUT : ",I8)') nkpt_
  write(*,*)
  stop
end if
if (nspin.ne.nspin_) then
  write(*,*)
  write(*,'("Error(getdmatdmft): differing nspinor")')
  write(*,'(" current    : ",I8)') nspin
  write(*,'(" DMATDMFT.OUT : ",I8)') nspin_
  write(*,*)
  stop
end if
!allocate the dmft density matrix
allocate(dmatkdmft(nstsv,nstsv,nkpt))
!zero density matrix
dmatkdmft(:,:,:)=0.d0
!set the diagonal to the KS/natural orbital occupancy values
do ist=1,nstsv
  !new occsv will be scaled by occmax after diagonalisation (in gwdmatk.f90)
  dmatkdmft(ist,ist,:)=occsv(ist,:)/occmax
enddo
!read in the density matrix
do ik=1,nkpt
  do ispn=1,nspin
    read(32,*) ik_,ispn_,minst,maxst
    allocate(rd(minst:maxst),id(minst:maxst))
    do ist=minst,maxst
      read(32,*) (rd(jst),id(jst),jst=minst,maxst)
      dmatkdmft(ist,minst:maxst,ik)=CMPLX(rd,id)
    enddo
    deallocate(rd,id)
  enddo
enddo 
!calculate the dmft temperature
tempk=1/(beta*kboltz)
if (mp_mpi) then
  write(*,*) 'DMFT temperature used: ',tempk,' K'
  write(*,*) 'DMFT chemical potential used: ',mu,' Hartree'
endif
!reads the hubbard energy minus the correction energy term in Hatree
read(32,*) engydmft    

close(32)
return
end subroutine 

