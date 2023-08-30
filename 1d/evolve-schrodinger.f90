program Gross_Pitaevskii_1d
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

  !integer :: i ! Automatic creation of array
  
  ! These are input parameters used for readability
  integer :: n, nf
  real(dl) :: lSize
  
  ! Move into time stepper
  integer :: n_out_steps
  real(dl) :: dt

  n = 512
  lSize = 80._dl

  nf = 3 ! Adjust based on PML choice
  call setup_simulation( nf, n, lSize/n, fld, tcur)
  call initialise_wavepacket(fld, 0._dl, 0._dl, 2._dl)
  !call initialise_coherent_state(fld, 0._dl, 0._dl)
  !call set_time_steps_decay()
  time_stepper%dt = twopi/64./8.
  time_stepper%out_size = 8
  time_stepper%n_out_steps = 1000
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
  subroutine initialise_coherent_state(fld, x0, p0)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: x0, p0

    fld = 0._dl
    fld(:,1,1) = exp(-0.5_dl*(xVals-x0)**2)*cos( p0*(xVals-0.5*x0) ) / (0.5_dl*twopi)**0.25
    fld(:,2,1) = exp(-0.5_dl*(xVals-x0)**2)*sin( p0*(xVals-0.5*x0) ) / (0.5_dl*twopi)**0.25
  end subroutine initialise_coherent_state 
  
  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup before starting the program
  subroutine setup_simulation(nf, nl, dx, fld, tcur)
    integer, intent(in) :: nf,nl
    real(dl), intent(in) :: dx
    real(dl), dimension(:,:,:), pointer :: fld
    real(dl), pointer :: tcur
    
    call create_lattice(nf,nl,dx)
    call init_integrator(nVar)

    fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
    tcur => yvec(2*nLat*nFld+1)
  end subroutine setup_simulation
  
  ! Convert this to take a TimeStepper.  Or better, a Model object
  subroutine time_evolve(dt,tend,dtout)
    real(dl), intent(in) :: dt
    real(dl), intent(in) :: tend, dtout
    integer :: outsize, nums
    integer :: i,j

    if (dt > dx**2/twopi) print*,"Warning, potentially undersampling the Nyquist mode"

    outsize = floor(dtout/dt)
    nums = floor(tend/dt)/outsize
    
    call output_fields_binary(fld)
    call output_log_file(fld, 0., dt)

    print*,"dt = ",dt," dt_out = ",dt*outsize, " num_out = ",nums
    
    tcur = 0._dl
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       call output_log_file(fld, dt*i*outsize, dt)
       call output_fields_binary(fld)
    enddo
  end subroutine time_evolve
  
end program Gross_Pitaevskii_1d
