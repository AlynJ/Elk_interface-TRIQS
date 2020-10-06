
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
subroutine hmlaa(ias,ngp,apwalm,ld,h)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients (in,complex(ngkmax,apwordmax,lmmaxapw))
!   ld     : leading dimension of h (in,integer)
!   h      : Hamiltonian matrix (inout,complex(*))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: h(*)
! local variables
integer is,io,jo
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
real(8) t0
complex(8) z1
! automatic arrays
complex(8) x(ngp),y(ngp)
is=idxis(ias)
t0=0.5d0*rmt(is)**2
lm1=0
do l1=0,lmaxapw
  do m1=-l1,l1
    lm1=lm1+1
    do io=1,apword(l1,is)
      y(:)=0.d0
      lm3=0
      do l3=0,lmaxapw
        do m3=-l3,l3
          lm3=lm3+1
          do jo=1,apword(l3,is)
            z1=0.d0
            do l2=0,lmaxo
              if (mod(l1+l2+l3,2).eq.0) then
                do m2=-l2,l2
                  lm2=idxlm(l2,m2)
                  z1=z1+gntyry(lm2,lm3,lm1)*haa(lm2,jo,l3,io,l1,ias)
                end do
              end if
            end do
            if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
              call zaxpy(ngp,z1,apwalm(:,jo,lm3),1,y,1)
            end if
          end do
        end do
      end do
! kinetic surface contribution
      do jo=1,apword(l1,is)
        z1=t0*apwfr(nrmt(is),1,io,l1,ias)*apwdfr(jo,l1,ias)
        call zaxpy(ngp,z1,apwalm(:,jo,lm1),1,y,1)
      end do
      x(1:ngp)=conjg(apwalm(1:ngp,io,lm1))
      call zher2i(ngp,x,y,ld,h)
    end do
  end do
end do
return

contains

pure subroutine zher2i(n,x,y,ld,a)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: x(n),y(n)
integer, intent(in) :: ld
complex(8), intent(inout) :: a(*)
! local variables
integer i,j,k
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) z1
do j=1,n
  k=(j-1)*ld
  z1=y(j)
  if (abs(dble(z1)).gt.eps) then
    if (abs(aimag(z1)).gt.eps) then
! complex prefactor
      do i=1,j-1
        a(k+i)=a(k+i)+z1*x(i)
      end do
      a(k+j)=dble(a(k+j))+dble(z1*x(j))
    else
! real prefactor
      a1=dble(z1)
      do i=1,j-1
        a(k+i)=a(k+i)+a1*x(i)
      end do
      a(k+j)=dble(a(k+j))+a1*dble(x(j))
    end if
  else if (abs(aimag(z1)).gt.eps) then
! imaginary prefactor
    b1=aimag(z1)
    do i=1,j-1
      a(k+i)=a(k+i)+b1*cmplx(-aimag(x(i)),dble(x(i)),8)
    end do
    a(k+j)=dble(a(k+j))-b1*aimag(x(j))
  end if
end do
return
end subroutine

end subroutine
!EOC

