!Written by A. D. N. James (copied from bandstr.f90)

! !INTERFACE:
subroutine writeband
use modmain
implicit none
! local variables
integer ist,iv,ik
real(8) emin,emax,e

!write out the band energies in normal format 
emin=1.d5
emax=-1.d5
open(50,file='BAND.OUT',form='FORMATTED')
do ist=1,nstsv
  do ik=1,nkpt
    e=evalsv(ist,ik)-efermi
    write(50,'(2G18.10)') dpp1d(ik),e
    emin=min(emin,e)
    emax=max(emax,e)
  end do
  write(50,'("     ")')
end do
close(50)
write(*,*)
write(*,'("Info(writewanproj):")')
write(*,'(" band structure plot written to BAND.OUT")')
open(50,file='BANDLINES.OUT',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),emin
  write(50,'(2G18.10)') dvp1d(iv),emax
  write(50,'("     ")')
end do
close(50)
write(*,*)
write(*,'(" vertex location lines written to BANDLINES.OUT")')
return
end subroutine

