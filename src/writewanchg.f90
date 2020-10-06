!Written by A. D. N. James 
subroutine writewanchg(ld,dmat)
use modmain
implicit none
! input
integer, intent(in) :: ld
complex(8), intent(inout) :: dmat(ld,nspinor,ld,nspinor,natmtot)
! local variables
integer ispn,jspn,l,lm1,lm2,i,j,iorb
integer lmmax,ia,is,ias
complex(8) chg,z
character(50) exfmt

!output the Wannier charge to file
open(810,file='WANCHARGE.OUT',form='FORMATTED')
do iorb=1,norb
! output charge in irreducible form
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  i=l**2+1
  j=(l+1)**2
  chg=0.d0
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! write to file
    write(810,'("Wan charge (Trace) of projected orbital ",I2,", atom ",I3," ")') iorb, ia
    write(810,*) ''
    write(exfmt,'("(",I8,"(F16.8))")') lmmax
    do ispn=1,nspinor
      do jspn=1,nspinor
        if((ispn.ne.jspn).and.(.not.spcpl)) cycle
        z=0.d0
        write(810,'("Wan density matrix of spinor block ",I2,",",I2,":")') ispn, jspn
        do lm1=i,j
          write(810,exfmt) (dble(dmat(lm1,ispn,lm2,jspn,ias)),lm2=i,j)
          z=z+dmat(lm1,ispn,lm1,jspn,ias)
        end do
        chg=chg+z
        write(810,'("Wan charge of spinor block = ",G14.8,"")') dble(z)
        write(810,'("")')
      end do
    end do
  end do
  write(810,'("Total wan charge = ",G14.8,"")') dble(chg)
  write(810,'("")')
  write(810,'("")')
end do

close(810)

return
end subroutine

