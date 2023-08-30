!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define FOURIER T

module lattice
  use constants, only : dl, twopi
  use utils, only : newunit
#ifdef FOURIER
  use fftw3
#endif
  implicit none
  
  ! Properties of simulation grid
  real(dl) :: dx, dk, len
  integer :: nFld, nLat
  integer :: nVar
  real(dl), dimension(:), allocatable, target :: yvec
  real(dl), dimension(:), allocatable :: pot
  real(dl), dimension(:), allocatable :: xVals
  real(dl), dimension(:), allocatable :: damping
#ifdef FOURIER
  type(transformPair1D) :: tPair 
#endif

  real(dl) :: del, nu 
  real(dl) :: mu

  type Lattice_1d
     real(dl), dimension(:), allocatable :: pot, xVals
     real(dl), dimension(:), allocatable :: damping
     real(dl) :: dx, dk, len
     integer :: nFld, nLat
     type(transformPair1D) :: tPair
  end type Lattice_1d
  
contains
  
  !>@brief
  !> Allocate storage for the lattice and current time
  !> Setup FFT for computing derivatives.
  subroutine create_lattice(nf,nl,dx_)
    integer, intent(in) :: nf, nl
    real(dl), intent(in) :: dx_

    integer :: i, u
    
    nFld = nf; nLat = nl
    nVar = 2*nFld*nLat+1; allocate( yvec(1:nVar) )
    allocate( pot(1:nLat) )
    allocate( xVals(1:nLat) )
    allocate( damping(1:nLat) )
    
#ifdef FOURIER
    call initialize_transform_1d(tPair,nLat)
#endif
    dx = dx_; len=nLat*dx; dk = twopi/len

    do i=1,nLat
       xVals(i) = (i-1)*dx - 0.5*len
    enddo

    call initialise_potential()
    call initialise_pml(2.*twopi, 2)

    !Output the potential
    open(unit=newunit(u), file='potential.dat')
    do i=1,nLat
       write(u,*) xVals(i), pot(i), damping(i)
    enddo
    close(u)
  end subroutine create_lattice

  ! Allow for some input here
  subroutine initialise_potential()
    integer :: i
    real(dl) :: x0
    real(dl) :: xc
    
    pot(:) = 0._dl 
    !pot(:) = 0.5_dl*xVals(:)**2

    ! Hertzberg
    !x0 = 2.75
    !pot(:) = 0.5_dl*xVals(:)**2
    !pot(:) = pot(:) * (1.-0.5*(xVals(:)/x0)**2) / (1.+0.5*(xVals(:)/x0)**4) + 0.5_dl

    ! Truncated Drummond
    do i=1,nLat
       xc = xVals(i)
       if ( abs(xc) <= 0.5_dl*twopi) then
          pot(i) = cos(xc) + 1._dl + 0.5_dl*1.4**2*sin(xc)**2
       endif
    enddo

  end subroutine initialise_potential

  subroutine initialise_pml(wid, pow)
    real(dl), intent(in) :: wid
    integer, intent(in) :: pow
    integer :: i
    real(dl) :: xc
    
    damping(:) = 0._dl

    do i=1,nLat
       xc = abs(xVals(i))
       if (xc >= wid) then
          !damping(i) = 100.*(abs(xVals(i))-wid)**2/32.**2
          damping(i) = (xc-wid)**pow / (1. + (xc-wid)**pow)
       endif
    enddo
  end subroutine initialise_pml
  
end module lattice
