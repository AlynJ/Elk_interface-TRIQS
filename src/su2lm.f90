subroutine  su2lm(lmmax,l,rlm,ias,subulm,sublm)

use modmain
use modmpi
!this routine calculates the SU(2) submatrix for the l input.
implicit none
!input
integer, intent(in) :: lmmax
integer, intent(in) :: l
integer, intent(in) :: rlm
integer, intent(in) :: ias
complex(8), intent(out) :: subulm(lmmax,lmmax)
integer, intent(inout) :: sublm(lmmax)
!local
real(8), allocatable :: elm(:,:)
integer ll, ld, lm, lm1, lm2, lmax, m, is, ia
complex(8), allocatable :: ulm(:,:,:),z(:),zz(:)
character(50) fname

!Output projectors are in the Spherical harmonic basis, but the
!projectors
!Will have been constructed in the user specified irreducible basis
! rotation matrix to real spherical harmonics (cubic basis)

! indices for desired portion of the rotation matrix
if(l.eq.0) then
  lm1=1
  lm2=1
! sub SU(2) matrix index (for l=0)
  sublm(1)=1
else
! max lm index for (l-1) to define stating lm index
  lm1=idxlm(l-1,l-1)+1
! max lm index for (l) 
  lm2=idxlm(l,l)
! sub SU(2) matrix indices (for l!=0)
  sublm(1:rlm)=sublm(1:rlm)-lm1+1
endif
!for cubic basis
if (cubic) then
  call rotcubic(l,lmmax,subulm)
elseif (lmirep) then
  lmax=min(3,lmaxo)
  ld=(lmax+1)**2
  allocate(elm(ld,natmtot),ulm(ld,ld,natmtot))
! generate SU(2)
  call genlmirep(lmax,ld,elm,ulm)
! irreducible representations file 
  is=idxis(ias)
  ia=idxia(ias)
  if (mp_mpi) then
    write(fname,'("_S",I2.2,"_A",I4.4,"")') is, ia
    open(50,file='ELMIREP'//trim(fname)//trim(filext),form='FORMATTED')
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    do ll=0,lmaxdos
     do m=-ll,ll
       lm=idxlm(ll,m)
       write(50,'(" l = ",I2,", m = ",I2,", lm= ",I3," : ",G18.10)') ll,m, &
                                                             lm,elm(lm,ias)
      end do
    end do
  endif
! for sub SU(2) matrix 
  subulm(1:lmmax,1:lmmax)=ulm(lm1:lm2,lm1:lm2,ias)
  deallocate(elm,ulm)
  close(50)
else
  !set subulm to the identity
  do lm=1,lmmax
    subulm(lm,lm)=1.d0
  enddo
end if
return

contains

subroutine rotcubic(l,lmmax,cub)
implicit none
!arguments
integer, intent(in) :: l
integer, intent(in) :: lmmax
complex(8), intent(inout) :: cub(lmmax,lmmax) !assumed to be zero at entry of subroutine
!local
integer i
!define the cubic unitary rotation matrix
if(l.eq.0) then
 cub(1,1)=1.d0
elseif(l.eq.1) then
 cub(1,1)=cmplx(1.d0/sqrt(2.d0),0.d0,8); cub(1,3)=cmplx(0.d0,-1.d0/sqrt(2.d0),8)
 cub(2,1)=cmplx(0.d0,1.d0/sqrt(2.d0),8); cub(2,3)=cmplx(0.d0,1.d0/sqrt(2.d0),8)
 cub(3,2)=cmplx(1.d0,0.d0,8)
elseif(l.eq.2) then
 cub(1,3)=cmplx(1.d0,0.d0,8)
 cub(2,1)=cmplx(1.d0/sqrt(2.d0),0.d0,8); cub(2,5)=cmplx(1.d0/sqrt(2.d0),0.d0,8)
 cub(3,1)=cmplx(-1.d0/sqrt(2.d0),0.d0,8); cub(3,5)=cmplx(1.d0/sqrt(2.d0),0.d0,8)
 cub(4,2)=cmplx(1.d0/sqrt(2.d0),0.d0,8); cub(4,4)=cmplx(-1.d0/sqrt(2.d0),0.d0,8)
 cub(5,2)=cmplx(1.d0/sqrt(2.d0),0.d0,8); cub(5,4)=cmplx(1.d0/sqrt(2.d0),0.d0,8)
elseif(l.eq.3) then
 cub(1,1)=cmplx(1.d0/sqrt(2.d0),0.d0,8); cub(1,7)=cmplx(-1.d0/sqrt(2.d0),0.d0,8)
 cub(2,2)=cmplx(1.d0/sqrt(2.d0),0.d0,8); cub(2,6)=cmplx(1.d0/sqrt(2.d0),0.d0,8)
 cub(3,3)=cmplx(1.d0/sqrt(2.d0),0.d0,8); cub(3,5)=cmplx(-1.d0/sqrt(2.d0),0.d0,8)
 cub(4,4)=cmplx(1.d0,0.d0,8)
 cub(5,3)=cmplx(0.d0,1.d0/sqrt(2.d0),8); cub(5,5)=cmplx(0.d0,1.d0/sqrt(2.d0),8)
 cub(6,2)=cmplx(0.d0,1.d0/sqrt(2.d0),8); cub(6,6)=cmplx(0.d0,-1.d0/sqrt(2.d0),8)
 cub(7,1)=cmplx(0.d0,1.d0/sqrt(2.d0),8); cub(7,7)=cmplx(0.d0,1.d0/sqrt(2.d0),8)
else
 write(*,*) 'Can not rotate spherical harmonic basis to cubic basis for l>3.'
 write(*,*) 'Spherical harmonics will be used instead'
 do i=1,lmmax
   cub(i,i)=1.d0
 enddo
endif
return
end subroutine

end subroutine
!EOC
