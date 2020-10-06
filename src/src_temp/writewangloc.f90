!Written by A. D. N. James 
subroutine writewangloc(tspndg,ld,nproj,gloc)

use modmain
implicit none
!arguments
logical, intent(in) :: tspndg
integer, intent(in) :: ld
integer, intent(in) :: nproj
complex(8), intent(in) :: gloc(nwplot,ld,nspinor,ld,nspinor,norb,nproj)
!local
integer ik,ist,jst,iw,ispn,jspn,lm1,lm2,iorb
integer nthd,is,ia,ias,iv,i,j,l,lmmax
real(8) eig,t1,dw,sf,emax,emin
complex(8), allocatable :: z1(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
real(8), allocatable :: wr(:)
character(50) exfmt,fname,fext

allocate(wr(nwplot))
dw=(wplot(2)-wplot(1))/dble(nwplot)
do iw=1,nwplot
  wr(iw)=dw*dble(iw-1)+wplot(1)
end do
!output the local Wannier spectral function
do iorb=1,norb
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do lm1=1,lmmax
      write(fname,'("GLOCPROJ_L",I2.2,"_S",I2.2,"_A",I4.4,"_",I2.2,"")') l,is,ia,lm1
      open(801,file=trim(fname)//trim(filext),form='FORMATTED')
      write(exfmt,'("(",I8,"(F16.8))")') lmmax+1
        do ispn=1,nspinor
          do jspn=1,nspinor
            if (tspndg.and.(ispn.ne.jspn)) cycle
            do iw=1,nwplot
            !write the spectral function
              write(801,exfmt)wr(iw), &
              (-aimag(gloc(iw,lm1,ispn,lm2,jspn,iorb,ia))/pi,lm2=1,lmmax)
          end do
          write(801,'("")')
        end do
      end do
    end do
    close(801)
  end do
end do
deallocate(wr)
return
end subroutine

