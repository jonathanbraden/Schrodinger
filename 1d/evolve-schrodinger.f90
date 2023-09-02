program Schrodinger_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use utils, only : newunit
  use eom
  use TimeStepping
  use output
  use integrator

  implicit none

  real(dl), dimension(:,:,:), pointer :: fld
  real(dl), pointer :: tcur
  type(TimeStepper) :: time_stepper
  
  ! These are input parameters used for readability
  integer :: n, nf
  real(dl) :: lSize
  
  ! Move into time stepper
  integer :: n_out_steps
  real(dl) :: dt

  n = 256
  lSize = 80._dl ! Adjust to x0 choice ...

  ! Note to self:
  !  Need to scale width of Gaussian to sigma^2 = 1/sqrt{m_FV^2}
  !  Oops, is it this way or the other way?
  
  nf = 3 ! Adjust based on PML choice
  call setup_simulation( nf, n, lSize/n, fld, tcur)
  !call initialise_wavepacket(fld, -20._dl, 5._dl, 1._dl)
  call initialise_coherent_state(fld, 0._dl, 0._dl, m2=1._dl)
  
  time_stepper%dt = twopi/64./8.
  time_stepper%out_size = 64
  time_stepper%n_out_steps = 8*100
  time_stepper%tcur = 0._dl
  call print_time_stepper(time_stepper)

  call time_evolve_stepper(fld, time_stepper, verbose_=.false.)
    
contains
  
  subroutine initialise_gaussian(fld, x0, sig)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl) :: x0, sig

    fld = 0._dl

    fld(:,1,1) = exp(-0.5_dl*(xVals-x0)**2/sig**2) / (0.5_dl*twopi*sig**2)**0.25
  end subroutine initialise_gaussian

  subroutine initialise_wavepacket(fld, x0, p0, sig)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: x0, p0, sig
    
    ! Fix this nonlocality with xVals
    ! Check ordering
    fld(:,2,:) = 0._dl
    fld(:,1,1) = exp(-0.5_dl*(xVals-x0)**2/sig**2)*cos( p0*xVals ) / (0.5_dl*twopi*sig**2)**0.25
    fld(:,2,1) = exp(-0.5_dl*(xVals-x0)**2/sig**2)*sin( p0*xVals ) / (0.5_dl*twopi*sig**2)**0.25
  end subroutine initialise_wavepacket

  ! Check which indexing is real/imaginary, and which is the field
  subroutine initialise_coherent_state(fld, x0, p0, m2)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: x0, p0
    real(dl), intent(in), optional :: m2

    real(dl) :: sig2

    sig2 = 1._dl
    if (present(m2)) then; sig2 = sqrt(m2); endif
       
    fld = 0._dl
    fld(:,1,1) = exp(-0.5_dl*(xVals-x0)**2/sig2)*cos( p0*(xVals-0.5*x0) ) / (0.5_dl*twopi*sig2)**0.25
    fld(:,2,1) = exp(-0.5_dl*(xVals-x0)**2/sig2)*sin( p0*(xVals-0.5*x0) ) / (0.5_dl*twopi*sig2)**0.25
  end subroutine initialise_coherent_state 
  
  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup before starting the program
  subroutine setup_simulation(nf, nl, dx, fld, tcur)
    integer, intent(in) :: nf,nl
    real(dl), intent(in) :: dx
    real(dl), dimension(:,:,:), pointer :: fld
    real(dl), pointer :: tcur

    real(dl) :: x0, pml_loc
    
    call create_lattice(nf,nl,dx)
    call init_integrator(nVar)
    x0 = 2.75_dl
    pml_loc = 4.*x0
    !pml_loc = 2.*twopi*x0
    call initialise_model(x0, pml_loc)

    !call initialise_model(x0, 2.*twopi*x0)  # BEC
    
    fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
    tcur => yvec(2*nLat*nFld+1)
  end subroutine setup_simulation
  
end program Schrodinger_1d
