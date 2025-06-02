!> @brief Deallocate and reset variables used for IC extrusion.
!>
!> @details
!>   Deallocates:
!>     - stored_values(:,:,:).
!>     - x_coords(:).
!>     - y_coords(:).
!>
!>   Resets:
!>     - files_loaded => .false.
!>     - index_x, index_y, global_offset_x => 0.
!>
!>   (Note: The corresponding loader macro reads
!>    `prim.<variable>.00.<timestep>.dat` files—pattern controlled by `zeros_default`—
!>    from `examples/Case_File/D`.)

#:def HardcodedDimensionsExtrusion()
    integer :: xRows, yRows, nRows, iix, iiy
    integer :: f, iter, ios, ios2, unit, unit2, idx, idy, index_x, index_y, jump, line_count, ycount
    real(wp) :: x_len, x_step, y_len, y_step
    real(wp) :: dummy_x, dummy_y, dummy_z, x0, y0
    integer :: global_offset_x, global_offset_y           ! MPI subdomain offset
    real(wp) :: delta_x, delta_y
    character(len=100), dimension(sys_size) :: fileNames ! Arrays to store all data from files
    character(len=200) :: errmsg
    real(wp), allocatable :: stored_values(:, :, :)
    real(wp), allocatable :: x_coords(:), y_coords(:)
    logical :: files_loaded = .false.
    real(wp) :: domain_xstart, domain_xend, domain_ystart, domain_yend
    character(len=*), parameter :: init_dir = "/home/YourDirectory" ! For example /home/MFC/examples/1D_Shock/D/
    character(len=20) :: file_num_str     ! For storing the file number as a string
    character(len=20) :: zeros_part       ! For the trailing zeros part
    character(len=6), parameter :: zeros_default = "00000"  ! Default zeros (can be changed)
#:enddef

#:def HardcodedDellacation()
    if (allocated(stored_values)) then
        deallocate (stored_values)
        deallocate (x_coords)
    end if

    if (allocated(y_coords)) then
        deallocate (y_coords)
    end if
#:enddef
