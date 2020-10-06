!Written by A. D. N. James 
subroutine writewancharge(ld,dmat)
use modmain
implicit none
! input
integer, intent(in) :: ld
complex(8), intent(in) ::dmat(ld,nspinor,ld,nspinor,natmtot)
! local variables
integer ist,jst,ispn,jspn,l,lm1,lm2,ik,i,j,iorb
integer ind1,ind2,lmmax,ia,ja,is,ias,jas,mst
complex(8), allocatable :: z1(:,:), z2(:,:)
complex(8), allocatable :: z3(:,:)
complex(8) chg,z
complex(8) su2(2,2),b(2,2),c(2,2)
real(8) f
character(50) exfmt

!output the Wannier charge to file
open(810,file='WANCHARGE.OUT',form='FORMATTED')
do iorb=1,norb
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  i=l**2+1
  j=(l+1)**2
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(810,'("Charge of Wannier projected atom ",I3," of species ",I3,"")') ia,is
    write(810,*) ''
    write(exfmt,'("(",I8,"(F16.8))")') lmmax
    chg=0.d0
    do ispn=1,nspinor
      do jspn=1,nspinor
        if((ispn.ne.jspn).and.(.not.spcpl)) cycle
        z=0.d0
        do lm1=i,j
          write(810,exfmt) (dble(dmat(lm1,ispn,lm2,jspn,ias)),lm2=i,j)
          z=z+dmat(lm1,ispn,lm1,jspn,ias)
        end do
        chg=chg+z
        write(810,'("")')
        write(810,'("Charge of projected orbitals of spinor",I2,",",I2," = ",G14.8,"")') &
                                                                     ispn, jspn, dble(z)
        write(810,'("")')
      end do
    end do
    write(810,'("Total Charge (Trace) of the lm orbitals = ",G14.8,"")') dble(chg)
    write(810,*) ''
  end do
end do
close(810)
return
end subroutine

