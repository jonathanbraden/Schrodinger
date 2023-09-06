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

    ! Temporary stuff to remove and pass in as parameters
    real(dl) :: x0
    
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
  end subroutine create_lattice

  subroutine initialise_model(x0, pml_loc)
    real(dl), intent(in) :: x0
    real(dl), intent(in) :: pml_loc
    integer :: u, i

    integer :: type

    type = 2

    call initialise_potential(x0, type)
    call initialise_pml(pml_loc, 2)

    open(unit=newunit(u), file='potential.dat')
    write(u,*) "# \mu x    V(x)   Gamma(x)"
    do i=1,nLat
       write(u,*) xVals(i), pot(i), damping(i)
    enddo
    close(u)
  end subroutine initialise_model
  
  ! Allow for some input here
  subroutine initialise_potential(x0, type)
    real(dl), intent(in) :: x0
    integer, intent(in) :: type
    
    integer :: i
    real(dl) :: xc, pot_max, x_max

    select case (type)

    case (1)
       pot(:) = 0.5_dl*xVals(:)**2

    case (2)
       x_max = sqrt(sqrt(3._dl)-1._dl)
       pot_max = 0.5_dl*( x_max**2*(1._dl-0.5_dl*x_max**2)/(1._dl+0.5_dl*x_max**4) + 1._dl )*x0**2

       pot(:) = 0.5_dl*xVals(:)**2
       pot(:) = pot(:) * (1._dl-0.5_dl*(xVals(:)/x0)**2) / (1._dl+0.5_dl*(xVals(:)/x0)**4) + 0.5_dl*x0**2

    case (3)
       x_max = 0._dl    ! Fill in
       pot_max = 0._dl  ! Fill in
       pot(:) = cos(xVals(:)/x0) + 1._dl + 0.5_dl*1.4**2*sin(xVals(:)/x0)**2
       pot(:) = x0**2 * pot(:)
       where (abs(xVals)>0.5_dl*twopi*x0) pot = 0._dl
       
    case default
       pot(:) = 0._dl

    end select
       
    ! Hertzberg
    !x_m = sqrt(sqrt(3._dl)-1._dl)
    !pot_m = (0.5_dl*x_m**2*(1.-0.5*x_m**2)/(1.+0.5*x_m**4) + 0.5_dl)*x0**2
    
    !pot(:) = 0.5_dl*xVals(:)**2
    !pot(:) = pot(:) * (1.-0.5*(xVals(:)/x0)**2) / (1.+0.5*(xVals(:)/x0)**4) + 0.5_dl*x0**2

    ! Modify at potential maximum
    !where (xVals(:) < -x_m*x0) pot = pot_m
    !where (xVals < -x_m*x0) pot = pot_m - (xVals+x_m*x0)**3/x0**3
    !where (pot > 16.) pot = 16.
    
    ! Extend quadratic all the way
    !where (xVals < 0.) pot = 0.5*xVals(:)**2 + 0.5_dl*x0**2

    call truncate_potential(16.)
  end subroutine initialise_potential

  subroutine extend_maximum(x_max, pot_max, x_scl, order)
    real(dl), intent(in) :: x_max, pot_max, x_scl
    integer, intent(in) :: order
    
    select case (order)

    case (0)
       where (xVals<x_max) pot = pot_max

    case (1)
       where (xVals<x_max) pot = pot_max - (xVals+x_max)/x_scl 

    case (3)
       where (xVals<x_max) pot = pot_max - ((xVals+x_max)/x_scl)**3

    case default
       print*,"Warning, only constant, linear and cubic deformations allowed in extend_maximum."
       print*,"Defaulting to constant"
       where (xVals<x_max) pot = pot_max
 
    end select    
  end subroutine extend_maximum
  
  ! Extend potential as quadratic to left of false vacuum
  subroutine extend_quadratic(x_fv, pot_fv)
    real(dl), intent(in) :: x_fv, pot_fv
    where (xVals < x_fv) pot(:) = 0.5*x_fv**2 + pot_fv
  end subroutine extend_quadratic

  ! Add option for a smooth truncation
  subroutine truncate_potential(thresh)
    real(dl), intent(in) :: thresh
    where (pot > thresh) pot = thresh
  end subroutine truncate_potential
  
  subroutine initialise_pml(wid, pow)
    real(dl), intent(in) :: wid
    integer, intent(in) :: pow
    integer :: i
    real(dl) :: xc
    
    damping(:) = 0._dl
    where (abs(xVals)>=wid) damping = (abs(xVals)-wid)**pow / (1._dl+(abs(xVals)-wid)**pow)
    !where (xVals>=wid) damping = (abs(xVals)-wid)**pow / (1._dl+(abs(xVals)-wid)**pow)
    
  end subroutine initialise_pml
  
end module lattice
