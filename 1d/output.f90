module Output
  use constants, only : dl
  use utils, only : newunit
  use eom ! needed for model parameters, etc.  Refactor to remove
  implicit none

contains
  
  ! Debug this so that all the parameters are passed in as arguments
  subroutine write_header(oFile, dt)
    integer, intent(in) :: oFile
    real(dl), intent(in) :: dt

    write(oFile,*) "# Schrodinger Simulation Lattice Parameters :"
    write(oFile,*) "# n_lattice = ",nLat," dx = ",dx
    write(oFile,*) "# Time Stepping Parameters :"
    write(oFile,*) "# dt = ",dt
    write(oFile,*) "# Model Parameters :"
    write(oFile,*) "# Not yet implemented"
    write(oFile,*) "#============================"
    write(oFile,*) "# m_Qt   mt   <|psi|^2>  <E>   <E>_by parts"
  end subroutine write_header
  
  subroutine output_log_file(fld, t, dt, fName)
    real(dl), dimension(:,:,:), intent(in) :: fld
    real(dl), intent(in) :: t, dt
    character(80), intent(in), optional :: fName

    character(80) :: fn
    integer, save :: oFile
    logical :: o
    real(dl) :: en, en_lap, prob

    fn = 'log.out'
    if (present(fName)) fn = trim(fName)

    inquire(file=trim(fn), opened=o)
    if (.not.o) then
       open(unit=newunit(oFile), file=trim(fn))
       call write_header(oFile, dt)
    endif

    en = compute_total_energy(fld)
    en_lap = compute_total_energy_laplacian(fld)
    prob = compute_probability(fld)
    write(oFile,*) t, 0._dl, prob, en, en_lap
    
  end subroutine output_log_file
  
  subroutine output_fields_binary(fld, fName)
    real(dl), dimension(:,:,:), intent(in) :: fld
    character(80), intent(in), optional :: fName

    logical :: o
    character(80) :: fn
    integer, save :: oFile

    fn = 'fields.bin'
    if (present(fName)) fn = trim(fName)

    inquire(file=trim(fn), opened=o)
    if (.not.o) open(unit=newunit(oFile), file=fn, access='stream', status='replace')

    write(oFile) fld
  end subroutine output_fields_binary
    
end module Output
