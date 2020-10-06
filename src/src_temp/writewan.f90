!Written by A. D. N. James

subroutine writewan
!This routine generates the Wannier projectors of all user input orbitals
use modmain
use moddftu

implicit none
!local variables
integer iorb,lm,ia,is,l,ja,ias,jas,lmmax,lmmax2,rlm
integer ld,isym,lspl,eqva,prv_atms,ik,iv,ist,nproj,nst
logical, allocatable :: done(:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: subulm(:,:,:,:)
integer, allocatable :: sublm(:,:,:)
complex(8), allocatable :: wanprj(:,:,:,:,:,:)
complex(8), allocatable :: gloc(:,:,:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
integer, allocatable :: idx(:,:), projst(:)
real(8), allocatable :: evalfv(:,:)

if(norb.le.0) then
  write(*,*) 'No projected orbital specified. Stopping.'
  stop
endif
!initalise the universal variables 
call init0
call init1
call init2
! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
!generate radial wavefunctions
call genapwfr
!generate LO radial wavefuntions
call genlofr
!calculate the radial overlap integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! generate the spin-orbit coupling radial functions
call gensocfr
!if outputting projectors into irreducible basis
if(cubic) lmirep=.true.
!if generating wannier bandstructure
if((task.eq.806).or.(task.eq.807).or.(task.eq.811)) then
! Output the energy eigenvalues and vectors to a new extension
  if(task.eq.807) then
    filext='_FS.OUT'
  elseif(task.eq.811) then
    filext='_WANBAND.OUT'
  endif
! generate desired eigenvalues and states and save them to
! binary files
  call genevfsv
! output energies and wkpt for non band structure tasks
  if(task.ne.811) then
    call writeeval
    call writekpts
    if(task.eq.806) call occupy
  else
    call writeband
  endif 
endif
!Calculate the desired projectors including equivalent atoms
!ld is max(2*l+1) of all projector inputs
ld=0
nproj=0
! rotation matrix to convert projector basis (initalised to zero)
do iorb=1,norb
  is = orb(iorb,1)
  l = orb(iorb,2)
  if(is.le.0) then
    write(*,*) '(writewanproj) incorrect species for orbital. Stopping.'
    stop
  endif
  if(l.lt.0) then
    write(*,*) '(writewanproj) incorrect l for orbital. Stopping.'
    stop
  endif
  ld=max(ld,2*l+1)
  nproj=max(nproj,natoms(is))
enddo
allocate(subulm(ld,ld,norb,nproj))
allocate(sublm(ld,norb,nproj))
subulm(:,:,:,:)=0.d0
sublm(:,:,:)=0.d0
lmmax=0
do iorb=1,norb
  is = orb(iorb,1)
  l = orb(iorb,2)
  lmmax=2*l+1
  do ia=1,natoms(is)
    rlm=orb(iorb,3)
    if((task.eq.805).or.(task.eq.806)) call writewansym(is,ia,l)
    ias=idxas(ia,is)
    sublm(1:rlm,iorb,ia)=rorblm(iorb,1:rlm)
    call su2lm(lmmax,l,rlm,ias,subulm(:,:,iorb,ia),sublm(:,iorb,ia))
    write(*,*) rlm,sublm(:,iorb,ia)
!check that the correct lm values have been inputted
    do lm=1,rlm
      if((sublm(lm,iorb,ia).le.0).or.(sublm(lm,iorb,ia).gt.lmmax)) then
        write(*,*) '(wannierproj) incorrect lm values input. Stopping'
        write(*,*) rlm,sublm(lm,iorb,ia)
        stop
      endif
    enddo
  enddo
enddo
!allocate the band indices in the correlated window
allocate(idx(nstsv,nkpt))
idx(:,:)=0
!allocate the number of band indices in the correlated window for each ik
allocate(projst(nkpt))
projst(:)=0
!Output energy window
write(*,'("Energy window:")') 
write(*,*) emincor, emaxcor
write(*,'("")') 
!determine the states in the band window
do ik=1,nkpt
!read the second variational energy eigenvalues from EVALSV.OUT
  if(task.eq.805) then
    call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
!put occupations for denser k-mesh projector calc
  elseif(task.eq.806) then
    call putoccsv(filext,ik,occsv)
  endif
! count and index states at k in energy window
  nst=0
!k index for reading in eigenvectors in genwfsvpwan
  do ist=1,nstsv
    if(((evalsv(ist,ik)-efermi).ge.emincor).and. &
       ((evalsv(ist,ik)-efermi).le.emaxcor)) then
      nst=nst+1
      idx(nst,ik)=ist
    end if
  end do
  projst(ik)=nst
enddo
!nst is now the maximum number of band indices
nst=maxval(projst(:))
!allocate the true k-dependent wannier projectors
allocate(wanprj(ld,nst,nspinor,norb,nproj,nkpt))
wanprj(:,:,:,:,:,:)=0.d0
! call wannier projector routines 
call wanproj(ld,nproj,nst,idx,projst,sublm,subulm,wanprj)
!Calculate local variables
if((task.ne.807).or.(task.ne.811)) then
!allocate the local density matrix
  lmmax2=(lmaxdm+1)**2
  allocate(dmat(lmmax2,nspinor,lmmax2,nspinor,natmtot))
  dmat(:,:,:,:,:)=0.d0
!calculate and output the wannier charge
  call wancharge(.false.,.false.,.false.,nproj,projst, &
               lmmax,lmmax2,nst,idx,subulm,wanprj,zzero,dmat)
  call writewancharge(lmmax2,dmat)
!calculate and output the wannier Green's functions
  allocate(gloc(nwplot,ld,nspinor,lmmax,nspinor,nproj))
  gloc(:,:,:,:,:,:)=0.d0
  call wangloc(.false.,.false.,.false.,nproj,projst, &
              ld,nst,idx,subulm,wanprj,gloc)
  call writewangloc(.false.,ld,nproj,gloc)
  deallocate(dmat,gloc)
endif
!write the projector files
!write out the general info for the projectors
call writewanflags(lmmax,nproj,subulm)
!write out the projectors
call writewanproj(lmmax,nst,nproj,projst,idx,wanprj)
deallocate(subulm,sublm,idx,wanprj,projst)
return
end subroutine
!EOC
