!Written by A. D. N. James

subroutine wannorm(nst,ld,nproj,sublm,wantemp,wantrue)
!Orthonormalises the temporary wannier projectors to calculate the true Wannier
!projectors
use modmain

implicit none
!arguments
integer, intent(in) :: nst
integer, intent(in) :: ld
integer, intent(in) :: nproj
integer, intent(in) :: sublm(ld,norb,nproj)
complex(8), intent(in) :: wantemp(ld,nst,nspinor,norb,nproj)
complex(8), intent(inout) :: wantrue(ld,nst,nspinor,norb,nproj)
! local variables
integer ispn,nspn,ist,lmfac,jspn,kspn,nd,iorb
integer lmmax,rlm,i,j,ia,is,l,lm,lm1
complex(8), allocatable :: D(:,:)
complex(8), allocatable :: olptmp(:,:)
complex(8), allocatable :: z(:,:), z1(:,:), z2(:,:), z3(:,:)

!need to normalise spinors of coupled spin systems by coupled spin
!normalisation matrix
if (spcpl) then
!number of spinors to loop over - set to one for coupled spins
  nspn=1
!couple spin elements in orthonormalisation
  lmfac=2
else
  nspn=nspinor
  lmfac=1
endif
!The dimension size of the overlap matrix [O = P(P*)] containing
!the overlap of all temporary projectors
nd=0
do iorb=1,norb
  is=orb(iorb,1) 
  rlm=orb(iorb,3) 
  do ia=1,natoms(is)
    nd=nd+rlm*lmfac
  end do
enddo
allocate(z(nd,nst),z1(nd,nst))
allocate(z2(nd,nd),z3(nd,nst))
allocate(olptmp(nd,nd),D(nd,nd))
!Calculate the Overlap matrix
do ispn=1,nspn
  olptmp(:,:)=0.d0
  D(:,:)=0.d0
!put temporary matrix into a 2d array to construct O.
  j=1
  do iorb=1,norb
    is=orb(iorb,1) 
    l=orb(iorb,2)
    rlm=orb(iorb,3)
    lmmax=2*l+1
    do ia=1,natoms(is)
      do lm=1,rlm
        lm1=sublm(lm,iorb,ia)
        z(j,:)=wantemp(lm1,:,ispn,iorb,ia)
        j=j+1
      enddo
!put the second spinor in the overlap matrix for spin
!coupled calculations
      if(spcpl) then
        do lm=1,rlm
          lm1=sublm(lm,iorb,ia)
          z(j,:)=wantemp(lm1,:,2,iorb,ia)
          j=j+1
        enddo
      endif
    enddo
  enddo
  z1(:,:)=z(:,:)
!calculate O = P(P*) 
  call zgemm('N','C',nd,nd,nst,zone, &
         z,nd,z,nd,zzero,olptmp,nd)
!determine the inverse square route of O (which is D)
  !get Z and Z*D^(1/2)
  call invsqrtwan(nd,olptmp,D)
  !calculate A^(1/2)=Z*D^(1/2)*Z^H
  !calculate A - Z=olptmp, Z*D^(1/2)=D
  z2(:,:)=0.d0
  call zgemm('N','C',nd,nd,nd,zone, &
         olptmp,nd,D,nd,zzero,z2,nd)
!calculate P_true = O^{-1/2}*P
  call zgemm('N','N',nd,nst,nd,zone, &
         z2,nd,z1,nd,zzero,z3,nd) 
! put into spinor projector
  j=1
  do iorb=1,norb
    is=orb(iorb,1) 
    l=orb(iorb,2)
    rlm=orb(iorb,3)
    lmmax=2*l+1
    do ia=1,natoms(is)
      do lm=1,rlm
        lm1=sublm(lm,iorb,ia)
        wantrue(lm1,:,ispn,iorb,ia)=z3(j,:)
        j=j+1
      enddo
!put the assign the second projector spinor for spin
!coupled calculations
      if(spcpl) then
        do lm=1,rlm
          lm1=sublm(lm,iorb,ia)
          wantrue(lm1,:,2,iorb,ia)=z3(j,:)
          j=j+1
        enddo
      endif
    enddo
  enddo
end do
deallocate(z,z1,z2,z3,D,olptmp)
return

contains

subroutine invsqrtwan(lmmax,Z,D)

!This subroutine calculates the square root of a Hermitian matrix A using the decomposition:
!                      A=Z*D*Z^H 
! where D is a diagonal matrix of eigenvalues of A,
! Z is matrix of orthonormal eigenvectors of A,
! Z^H is its Hermitian conjugate.
! Then A^(1/2)=Z*D^(1/2)*Z^H.                                     
implicit none
!arguments
integer, intent(in) :: lmmax
complex(8), intent(inout) :: Z(lmmax,lmmax)
complex(8), intent(inout) :: D(lmmax,lmmax)
!local
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)
real(8) w(lmmax)
complex(8) eigval(lmmax)
integer info,i,j,m

!Performing the LDA decomposition of the overlap matrix
allocate(work(2*lmmax-1))
allocate(rwork(3*lmmax-2))
call zheev('V', 'U', lmmax, Z, lmmax, w, work,2*lmmax-1, rwork, info)
!A is now the eigenvectors Z, making sure that the decomposition worked
if(info.ne.0) then
  write(*,*)'wanproj: zheev has returned a non zero info value - ',info
  write(*,*)'         The inverse square matrix could not be calculated.'
  write(*,*)'         This interface will stop.'
  stop
endif
!making sure that the calculated eigenvalues are okay for calculating the 
!projectors. Warns the user that the resulting projectors may be wrong 
do i=1,lmmax
  if(w(i).lt.0d0) then
    write(*,*)'wanproj WARNING: Eigenvalue of overlap matrix (', i,') is negative.'
    write(*,*)'                 The projectors may be wrong. Try using a different energy window'
  endif
  if(w(i).lt.1d-12) then
    write(*,*)'wanproj WARNING: Eigenvalue of overlap matrix (', i,') is close to 0.'
    write(*,*)'                 The projectors may be wrong. Try using a different energy window'
  endif
  if(w(i).eq.0.d0) then
    write(*,*)'wanproj: Eigenvalue of overlap matrix (', i,') is equal to 0.'
    write(*,*)'         Can not calculate the inverse square matrix.'
    write(*,*)'         Stopping.'
    stop
  endif  
enddo
!Calculating Z*D^{-1/2}
D=0.d0
eigval(:)=cmplx(w(:),0.d0)
do i=1,lmmax
  do j=1,lmmax
    D(i,j)=Z(i,j)/sqrt(eigval(j))
  end do
end do
return
end subroutine

end subroutine
!EOC
