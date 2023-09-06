!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the variables and equations describing the GPE we are solving
!>
!> This module provides storage for the complex fields describing a collection
!> of Bose condensates, along with appropriate time evolution routines.
!> We assume the condensates obey the Gross-Pitaevskii equation
!>  \f[
!>    i\hbar \frac{\partial\psi_i}{\partial t} = -\frac{\hbar^2}{2m}\nabla^2\psi_i - \mu_i\psi_i + \sum_j g_{ij}\left|\psi_j\right|^2\psi_i - \sum_j\nu_{ij}psi_j
!>  \f]
!> with the coefficients \f$g_{ij},\mu_i\f$ and \f$\nu_{ij}\f$ viewed as model parameters
!> to be adjusted in the experimental setup.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FOURIER T
!#define DISCRETE T

#define R1 1:nLat
#define I1 (nLat+1):(2*nLat)
#define R2 (2*nLat+1):(3*nLat)
#define I2 (3*nLat+1):(4*nLat)
#define R3 (4*nLat+1):(5*nLat)
#define I3 (5*nLat+1):(6*nLat)

module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  use lattice
  implicit none
    
contains

  subroutine derivs(yc, yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    yp = 0._dl  ! Remove this later

    ! Convert these to chebyshev derivatives
    tPair%realSpace(:) = yc(R1)
    call laplacian_1d_wtype(tPair,dk)
    yp(I1) = yp(I1) + 0.5_dl*tPair%realSpace(:) - pot(:)*yc(R1)
    
    tPair%realSpace(:) = yc(I1)
    call laplacian_1d_wtype(tPair,dk)
    yp(R1) = yp(R1) - 0.5_dl*tPair%realSpace(:) + pot(:)*yc(I1)

    yp(I1) = yp(I1) - damping(:)*yc(I1)
    yp(R1) = yp(R1) - damping(:)*yc(R1)
    yp(I2) = yp(I2) - damping(:)*yc(I2)
    yp(R2) = yp(R2) - damping(:)*yc(R2)

    yp(I1) = yp(I1) - damping(:)*pot(:)*yc(R3)
    yp(R1) = yp(R1) + damping(:)*pot(:)*yc(I3)
    
    tPair%realSpace(:) = yc(R1)
    call derivative_1d_wtype(tPair, dk)
    yp(I2) = yp(I2) - damping(:)*tPair%realSpace(:)

    tPair%realSpace(:) = yc(I1)
    call derivative_1d_wtype(tPair, dk)
    yp(R2) = yp(R2) + damping(:)*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(R2)
    call derivative_1d_wtype(tPair, dk)
    yp(R1) = yp(R1) + 0.5_dl*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(I2)
    call derivative_1d_wtype(tPair, dk)
    yp(I1) = yp(I1) + 0.5_dl*tPair%realSpace(:)

    yp(R3) = yc(R1)
    yp(I3) = yc(I1)
    
  end subroutine derivs
  
  !>@brief
  !> Compute the derivatives for Schrodinger equation
  subroutine derivs_cap(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    yp = 0._dl  ! Remove this later

    tPair%realSpace(:) = yc(R1)
    call laplacian_1d_wtype(tPair,dk)
    yp(I1) = yp(I1) + 0.5_dl*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(I1)
    call laplacian_1d_wtype(tPair,dk)
    yp(R1) = yp(R1) - 0.5_dl*tPair%realSpace(:)
    
    yp(I1) = yp(I1) - pot(:)*yc(R1)
    yp(R1) = yp(R1) + pot(:)*yc(I1)

    ! For testing purposes, put in CAP potential here
    !yp(R1) = yp(R1) - damping(:)*yc(R1)
    !yp(I1) = yp(I1) - damping(:)*yc(I1)
  end subroutine derivs_cap

  subroutine derivs_b(yc, yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    yp = 0._dl  ! Remove this later

    ! Convert these to chebyshev derivatives
    tPair%realSpace(:) = yc(R1)
    call laplacian_1d_wtype(tPair,dk)
    yp(I1) = yp(I1) + 0.5_dl*tPair%realSpace(:) - pot(:)*yc(R1)
    
    tPair%realSpace(:) = yc(I1)
    call laplacian_1d_wtype(tPair,dk)
    yp(R1) = yp(R1) - 0.5_dl*tPair%realSpace(:) + pot(:)*yc(I1)
    
    ! Add PML part
    yp(I1) = yp(I1) - damping(:)*yc(I1)
    yp(R1) = yp(R1) - damping(:)*yc(R1)
    yp(I2) = yp(I2) - damping(:)*yc(I2)
    yp(R2) = yp(R2) - damping(:)*yc(R2)

    ! Bug check these
    tPair%realSpace(:) = yc(R1)
    call derivative_1d_wtype(tPair, dk)
    yp(I2) = yp(I2) - damping(:)*tPair%realSpace(:)

    tPair%realSpace(:) = yc(I1)
    call derivative_1d_wtype(tPair, dk)
    yp(R2) = yp(R2) + damping(:)*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(R2)
    call derivative_1d_wtype(tPair, dk)
    yp(R1) = yp(R1) + 0.5_dl*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(I2)
    call derivative_1d_wtype(tPair, dk)
    yp(I1) = yp(I1) + 0.5_dl*tPair%realSpace(:)
  end subroutine derivs_b
  
  subroutine compute_probability_current(current, fld)
    real(dl), dimension(1:nLat), intent(out) :: current
    real(dl), dimension(:,:,:), intent(in) :: fld

    current = 0._dl
  end subroutine compute_probability_current

  subroutine compute_current_divergence(dprob, fld)
    real(dl), dimension(1:nLat), intent(out) :: dprob
    real(dl), dimension(:,:,:), intent(in) :: fld
  end subroutine compute_current_divergence
  
  ! Need to debug
  real(dl) function compute_total_energy(fld) result(en)
    real(dl), dimension(:,:,:), intent(in) :: fld

    integer :: j

    en = 0._dl
    
    do j=1,2 
       tPair%realSpace(:) = fld(:,j,1)
       call gradsquared_1d_wtype(tPair,dk)
       en = en + 0.5_dl*sum(tPair%realSpace)
       en = en + sum(pot(:)*fld(:,j,1)**2)
    enddo    
    en = en * dx
  end function compute_total_energy

  real(dl) function compute_total_energy_laplacian(fld) result(en)
    real(dl), dimension(:,:,:), intent(in) :: fld

    integer :: j

    en = 0._dl
    do j=1,2
       tPair%realSpace(:) = fld(:,j,1)
       call laplacian_1d_wtype(tPair, dk)
       en = en - 0.5_dl*sum( fld(:,j,1)*tPair%realSpace(:) )
       en = en + sum(pot(:)*fld(:,j,1)**2)
    enddo
    en = en * dx
  end function compute_total_energy_laplacian
       
  ! Fix normalizaion here
  real(dl) function compute_probability(fld) result(rho)
    real(dl), dimension(:,:,:), intent(in) :: fld

    rho = sum(fld(:,:,1)**2) * dx
  end function compute_probability
  
end module eom
