!>
!! @file
!! @author Henry Le Berre <hberre3@gatech.edu>
!! @brief  Contains module m_model

#:include 'macros.fpp'

!> @brief Binary STL file reader and processor for immersed boundary geometry
module m_model

    use m_helper
    use m_mpi_proxy
    use m_derived_types

    use iso_c_binding, only: c_char, c_int32_t, c_int16_t, c_float

    implicit none

    private

    public :: f_model_read, s_model_write, s_model_free, f_model_is_inside, models, gpu_ntrs, &
              gpu_trs_v, gpu_trs_n, gpu_boundary_v, gpu_boundary_edge_count, &
              gpu_total_vertices, stl_bounding_boxes

    ! Subroutines for STL immersed boundaries
    public :: f_check_boundary, f_register_edge, f_model_is_inside_flat, &
              f_distance_normals_3D, f_distance_normals_2D, s_pack_model_for_gpu

    !! array of STL models that can be allocated and then used in IB marker and levelset compute
    type(t_model_array), allocatable, target :: models(:)
    !! GPU-friendly flat arrays for STL model data
    integer, allocatable :: gpu_ntrs(:)
    real(wp), allocatable, dimension(:, :, :, :) :: gpu_trs_v
    real(wp), allocatable, dimension(:, :, :) :: gpu_trs_n
    real(wp), allocatable, dimension(:, :, :, :) :: gpu_boundary_v
    integer, allocatable :: gpu_boundary_edge_count(:)
    integer, allocatable :: gpu_total_vertices(:)
    real(wp), allocatable :: stl_bounding_boxes(:, :, :)

contains

    !> This procedure reads a binary STL file.
    !! @param filepath Path to the STL file.
    !! @param model The binary of the STL file.
    impure subroutine s_read_stl_binary(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(out) :: model

        integer :: i, iunit, iostat

        character(kind=c_char, len=80) :: header
        integer(kind=c_int32_t) :: nTriangles

        real(kind=c_float) :: normal(3), v(3, 3), v_norm
        integer(kind=c_int16_t) :: attribute

        open (newunit=iunit, file=filepath, action='READ', &
              form='UNFORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open Binary STL file ", filepath

            call s_mpi_abort()
        end if

        read (iunit, iostat=iostat) header, nTriangles

        if (iostat /= 0) then
            print *, "Error: could not read header from Binary STL file ", filepath

            call s_mpi_abort()
        end if

        model%ntrs = nTriangles

        allocate (model%trs(model%ntrs))

        do i = 1, model%ntrs
            read (iunit) normal(:), v(1, :), v(2, :), v(3, :), attribute

            model%trs(i)%v = v
            model%trs(i)%n = normal
            v_norm = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
            if (v_norm > 0._wp) model%trs(i)%n = normal/v_norm
        end do

        close (iunit)

    end subroutine s_read_stl_binary

    !> This procedure reads an ASCII STL file.
    !! @param filepath Path to the STL file.
    !! @param model the STL file.
    impure subroutine s_read_stl_ascii(filepath, model)
        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(out) :: model

        integer :: i, j, iunit, iostat
        character(80) :: line, buffered_line
        logical :: is_buffered
        real(wp) :: normal(3), v_norm

        is_buffered = .false.

        open (newunit=iunit, file=filepath, action='READ', &
              form='FORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open ASCII STL file ", filepath
            call s_mpi_abort()
        end if

        model%ntrs = 0
        do
            if (is_buffered) then
                line = buffered_line
                is_buffered = .false.
            else
                if (.not. f_read_line(iunit, line)) exit
            end if

            if (line(1:6) == "facet ") then
                model%ntrs = model%ntrs + 1
            end if
        end do

        allocate (model%trs(model%ntrs))

        rewind (iunit)

        i = 1
        do
            if (is_buffered) then
                line = buffered_line
                is_buffered = .false.
            else
                if (.not. f_read_line(iunit, line)) exit
            end if

            if (line(1:5) == "solid") cycle
            if (line(1:8) == "endsolid") exit

            if (line(1:12) /= "facet normal") then
                print *, "Error: expected facet normal in STL file ", filepath
                call s_mpi_abort()
            end if

            call s_skip_ignored_lines(iunit, buffered_line, is_buffered)
            read (line(13:), *) normal
            v_norm = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
            if (v_norm > 0._wp) model%trs(i)%n = normal/v_norm

            call s_skip_ignored_lines(iunit, buffered_line, is_buffered)
            if (is_buffered) then
                line = buffered_line
                is_buffered = .false.
            else
                read (iunit, '(A)') line
            end if

            do j = 1, 3
                if (is_buffered) then
                    line = buffered_line
                    is_buffered = .false.
                else
                    if (.not. f_read_line(iunit, line)) exit
                end if

                if (line(1:6) /= "vertex") then
                    print *, "Error: expected vertex in STL file ", filepath
                    call s_mpi_abort()
                end if

                call s_skip_ignored_lines(iunit, buffered_line, is_buffered)
                read (line(7:), *) model%trs(i)%v(j, :)
            end do

            if (is_buffered) then
                line = buffered_line
                is_buffered = .false.
            else
                if (.not. f_read_line(iunit, line)) exit
            end if

            if (is_buffered) then
                line = buffered_line
                is_buffered = .false.
            else
                if (.not. f_read_line(iunit, line)) exit
            end if

            if (line(1:8) /= "endfacet") then
                print *, "Error: expected endfacet in STL file ", filepath
                call s_mpi_abort()
            end if

            i = i + 1
        end do
    end subroutine s_read_stl_ascii

    !> This procedure reads an STL file.
    !! @param filepath Path to the STL file.
    !! @param model the STL file.
    impure subroutine s_read_stl(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(out) :: model

        integer :: iunit, iostat

        character(80) :: line

        open (newunit=iunit, file=filepath, action='READ', &
              form='FORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath

            call s_mpi_abort()
        end if

        read (iunit, '(A)') line

        close (iunit)

        if (line(1:5) == "solid") then
            call s_read_stl_ascii(filepath, model)
        else
            call s_read_stl_binary(filepath, model)
        end if

    end subroutine s_read_stl

    !> This procedure reads an OBJ file.
    !! @param filepath Path to the odj file.
    !! @param model The obj file.
    impure subroutine s_read_obj(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(out) :: model

        integer :: i, j, k, l, iunit, iostat, nVertices

        real(wp), dimension(1:3), allocatable :: vertices(:, :)

        character(80) :: line

        open (newunit=iunit, file=filepath, action='READ', &
              form='FORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open model file ", filepath

            call s_mpi_abort()
        end if

        nVertices = 0
        model%ntrs = 0
        do
            if (.not. f_read_line(iunit, line)) exit

            select case (line(1:2))
            case ("v ")
                nVertices = nVertices + 1
            case ("f ")
                model%ntrs = model%ntrs + 1
            end select
        end do

        rewind (iunit)

        allocate (vertices(nVertices, 1:3))
        allocate (model%trs(model%ntrs))

        i = 1
        j = 1

        do
            if (.not. f_read_line(iunit, line)) exit

            select case (line(1:2))
            case ("g ")
            case ("vn")
            case ("vt")
            case ("l ")
            case ("v ")
                read (line(3:), *) vertices(i, :)
                i = i + 1
            case ("f ")
                read (line(3:), *) k, l, j
                model%trs(j)%v(1, :) = vertices(k, :)
                model%trs(j)%v(2, :) = vertices(l, :)
                model%trs(j)%v(3, :) = vertices(j, :)
                j = j + 1
            case default
                print *, "Error: unknown line type in OBJ file ", filepath
                print *, "Line: ", line

                call s_mpi_abort()
            end select
        end do

        deallocate (vertices)

        close (iunit)

    end subroutine s_read_obj

    !> This procedure reads a mesh from a file.
    !! @param filepath Path to the file to read.
    !! @return The model read from the file.
    impure function f_model_read(filepath) result(model)

        character(LEN=*), intent(in) :: filepath

        type(t_model) :: model

        select case (filepath(len(trim(filepath)) - 3:len(trim(filepath))))
        case (".stl")
            call s_read_stl(filepath, model)
        case (".obj")
            call s_read_obj(filepath, model)
        case default
            print *, "Error: unknown model file format for file ", filepath

            call s_mpi_abort()
        end select

    end function f_model_read

    !> This procedure writes a binary STL file.
    !! @param filepath Path to the STL file.
    !! @param model STL to write
    impure subroutine s_write_stl(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(in) :: model

        integer :: i, j, iunit, iostat

        character(kind=c_char, len=80), parameter :: header = "Model file written by MFC."
        integer(kind=c_int32_t) :: nTriangles
        real(wp) :: normal(3), v(3)
        integer(kind=c_int16_t) :: attribute

        open (newunit=iunit, file=filepath, action='WRITE', &
              form='UNFORMATTED', iostat=iostat, access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath

            call s_mpi_abort()
        end if

        nTriangles = model%ntrs
        write (iunit, iostat=iostat) header, nTriangles

        if (iostat /= 0) then
            print *, "Error: could not write header to STL file ", filepath

            call s_mpi_abort()
        end if

        do i = 1, model%ntrs
            normal = model%trs(i)%n
            write (iunit) normal

            do j = 1, 3
                v = model%trs(i)%v(j, :)
                write (iunit) v(:)
            end do

            attribute = 0
            write (iunit) attribute
        end do

        close (iunit)

    end subroutine s_write_stl

    !> This procedure writes an OBJ file.
    !! @param filepath Path to the obj file.
    !! @param model obj to write.
    impure subroutine s_write_obj(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(in) :: model

        integer :: iunit, iostat

        integer :: i, j

        open (newunit=iunit, file=filepath, action='WRITE', &
              form='FORMATTED', iostat=iostat, access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open OBJ file ", filepath

            call s_mpi_abort()
        end if

        write (iunit, '(A)') "# Model file written by MFC."

        do i = 1, model%ntrs
            do j = 1, 3
                write (iunit, '(A, " ", (f30.20), " ", (f30.20), " ", (f30.20))') &
                    "v", model%trs(i)%v(j, 1), model%trs(i)%v(j, 2), model%trs(i)%v(j, 3)
            end do

            write (iunit, '(A, " ", I0, " ", I0, " ", I0)') &
                "f", i*3 - 2, i*3 - 1, i*3
        end do

        close (iunit)

    end subroutine s_write_obj

    !> This procedure writes a binary STL file.
    !! @param filepath  Path to the file to write.
    !! @param model Model to write.
    impure subroutine s_model_write(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(in) :: model

        select case (filepath(len(trim(filepath)) - 3:len(trim(filepath))))
        case (".stl")
            call s_write_stl(filepath, model)
        case (".obj")
            call s_write_obj(filepath, model)
        case default
            print *, "Error: unknown model file format for file ", filepath

            call s_mpi_abort()
        end select

    end subroutine s_model_write

    !> This procedure frees the memory allocated for an STL mesh.
    subroutine s_model_free(model)

        type(t_model), intent(inout) :: model

        deallocate (model%trs)

    end subroutine s_model_free

    impure function f_read_line(iunit, line) result(bIsLine)

        integer, intent(in) :: iunit
        character(80), intent(out) :: line

        logical :: bIsLine
        integer :: iostat

        bIsLine = .true.

        do
            read (iunit, '(A)', iostat=iostat) line

            if (iostat < 0) then
                bIsLine = .false.
                exit
            end if

            line = adjustl(trim(line))

            if (len(trim(line)) == 0) cycle
            if (line(1:5) == "solid") cycle
            if (line(1:1) == "#") cycle

            exit
        end do

    end function f_read_line

    !> @brief Reads the next non-comment line from a model file, using a buffered look-ahead mechanism.
    impure subroutine s_skip_ignored_lines(iunit, buffered_line, is_buffered)
        integer, intent(in) :: iunit
        character(80), intent(inout) :: buffered_line
        logical, intent(inout) :: is_buffered

        character(80) :: line

        if (is_buffered) then
            line = buffered_line
            is_buffered = .false.
        else
            if (.not. f_read_line(iunit, line)) return
        end if

        buffered_line = line
        is_buffered = .true.
    end subroutine s_skip_ignored_lines

    !> This function is used to replace the fortran random number
    !! generator because the native generator is not compatible being called
    !! from GPU routines/functions
    function f_model_random_number(seed) result(rval)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(inout) :: seed
        real(wp) :: rval

        seed = ieor(seed, ishft(seed, 13))
        seed = ieor(seed, ishft(seed, -17))
        seed = ieor(seed, ishft(seed, 5))

        rval = abs(real(seed, wp))/real(huge(seed), wp)
    end function f_model_random_number

    !> This procedure, recursively, finds whether a point is inside an octree.
    !! @param model    Model to search in.
    !! @param point    Point to test.
    !! @param spacing  Space around the point to search in (grid spacing).
    !! @param spc      Number of samples per cell.
    !! @return True if the point is inside the octree, false otherwise.
    impure function f_model_is_inside(model, point, spacing, spc) result(fraction)

        ! $:GPU_ROUTINE(parallelism='[seq]')

        type(t_model), intent(in) :: model
        real(wp), dimension(1:3), intent(in) :: point
        real(wp), dimension(1:3), intent(in) :: spacing
        integer, intent(in) :: spc
        real(wp) :: phi, theta
        integer :: rand_seed

        real(wp) :: fraction

        type(t_ray) :: ray
        integer :: i, j, k, nInOrOut, nHits

        real(wp), dimension(1:spc, 1:3) :: ray_origins, ray_dirs

        rand_seed = int(point(1)*73856093_wp) + &
                    int(point(2)*19349663_wp) + &
                    int(point(3)*83492791_wp)
        if (rand_seed == 0) rand_seed = 1

        ! generate our random collection or rays
        do i = 1, spc
            do k = 1, 3
                ! random jitter in the origin helps us estimate volume fraction instead of only at the cell center
                ray_origins(i, k) = point(k) + (f_model_random_number(rand_seed) - 0.5_wp)*spacing(k)
                ! cast sample rays in all directions
                ray_dirs(i, k) = point(k) + f_model_random_number(rand_seed) - 0.5_wp
            end do
            ray_dirs(i, :) = ray_dirs(i, :)/sqrt(sum(ray_dirs(i, :)*ray_dirs(i, :)))
        end do

        ! ray trace
        nInOrOut = 0
        do i = 1, spc
            ray%o = ray_origins(i, :)
            ray%d = ray_dirs(i, :)

            nHits = 0
            do j = 1, model%ntrs
                ! count the number of triangles this ray intersects
                if (f_intersects_triangle(ray, model%trs(j))) then
                    nHits = nHits + 1
                end if
            end do

            ! if the ray hits an odd number of triangles on its way out, then
            ! it must be on the inside of the model
            nInOrOut = nInOrOut + mod(nHits, 2)
        end do

        fraction = real(nInOrOut)/real(spc)

    end function f_model_is_inside

    impure function f_model_is_inside_flat(ntrs, trs_v, trs_n, pid, point, spacing, spc) result(fraction)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: ntrs
        real(wp), dimension(:, :, :, :), intent(in) :: trs_v
        real(wp), dimension(:, :, :), intent(in) :: trs_n
        integer, intent(in) :: pid
        real(wp), dimension(1:3), intent(in) :: point
        real(wp), dimension(1:3), intent(in) :: spacing
        integer, intent(in) :: spc

        real(wp) :: fraction
        real(wp) :: origin(1:3), dir(1:3), dir_mag
        type(t_ray) :: ray
        type(t_triangle) :: tri
        integer :: i, j, k, nInOrOut, nHits
        integer :: rand_seed

        rand_seed = int(point(1)*73856093_wp) + &
                    int(point(2)*19349663_wp) + &
                    int(point(3)*83492791_wp)
        if (rand_seed == 0) rand_seed = 1

        ! generate our random collection of rays
        nInOrOut = 0
        do i = 1, spc
            ! Generate one ray at a time — no arrays needed
            do k = 1, 3
                origin(k) = point(k) + (f_model_random_number(rand_seed) - 0.5_wp)*spacing(k)
                dir(k) = point(k) + f_model_random_number(rand_seed) - 0.5_wp
            end do
            dir_mag = sqrt(dir(1)*dir(1) + dir(2)*dir(2) + dir(3)*dir(3))
            dir(:) = dir(:)/dir_mag

            ray%o = origin
            ray%d = dir

            nHits = 0
            do j = 1, ntrs
                tri%v(:, :) = trs_v(:, :, j, pid)
                tri%n(:) = trs_n(:, j, pid)
                if (f_intersects_triangle(ray, tri)) then
                    nHits = nHits + 1
                end if
            end do
            nInOrOut = nInOrOut + mod(nHits, 2)
        end do

        fraction = real(nInOrOut)/real(spc)
    end function f_model_is_inside_flat

    ! From https://www.scratchapixel.com/lessons/3e-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
    !> This procedure checks if a ray intersects a triangle.
    !! @param ray      Ray.
    !! @param triangle Triangle.
    !! @return         True if the ray intersects the triangle, false otherwise.
    elemental function f_intersects_triangle(ray, triangle) result(intersects)

        $:GPU_ROUTINE(parallelism='[seq]')

        type(t_ray), intent(in) :: ray
        type(t_triangle), intent(in) :: triangle

        logical :: intersects

        real(wp) :: N(3), P(3), C(3), edge(3), vp(3)
        real(wp) :: area2, d, t, NdotRayDirection

        intersects = .false.

        N = triangle%n
        area2 = sqrt(sum(N(:)*N(:)))

        NdotRayDirection = sum(N(:)*ray%d(:))

        if (abs(NdotRayDirection) < 0.0000001_wp) then
            return
        end if

        d = -sum(N(:)*triangle%v(1, :))
        t = -(sum(N(:)*ray%o(:)) + d)/NdotRayDirection

        if (t < 0) then
            return
        end if

        P = ray%o + t*ray%d

        edge = triangle%v(2, :) - triangle%v(1, :)
        vp = P - triangle%v(1, :)
        C = f_cross(edge, vp)
        if (sum(N(:)*C(:)) < 0) then
            return
        end if

        edge = triangle%v(3, :) - triangle%v(2, :)
        vp = P - triangle%v(2, :)
        C = f_cross(edge, vp)
        if (sum(N(:)*C(:)) < 0) then
            return
        end if

        edge = triangle%v(1, :) - triangle%v(3, :)
        vp = P - triangle%v(3, :)
        C = f_cross(edge, vp)
        if (sum(N(:)*C(:)) < 0) then
            return
        end if

        intersects = .true.

    end function f_intersects_triangle

    !> This procedure checks and labels edges shared by two or more triangles facets of the 2D STL model.
    !! @param model                      Model to search in.
    !! @param boundary_vertex_count      Output total boundary vertex count
    subroutine f_check_boundary(model, boundary_v, boundary_vertex_count, boundary_edge_count)

        type(t_model), intent(in) :: model
        real(wp), allocatable, intent(out), dimension(:, :, :) :: boundary_v !< Output boundary vertices/normals
        integer, intent(out) :: boundary_vertex_count, boundary_edge_count !< Output boundary vertex/edge count

        integer :: i, j !< Model index iterator
        integer :: edge_count, edge_index, store_index !< Boundary edge index iterator
        real(wp), dimension(1:2, 1:2) :: edge !< Edge end points buffer
        real(wp), dimension(1:2) :: boundary_edge !< Boundary edge end points buffer
        real(wp), dimension(1:(3*model%ntrs), 1:2, 1:2) :: temp_boundary_v !< Temporary boundary vertex buffer
        integer, dimension(1:(3*model%ntrs)) :: edge_occurrence !< The manifoldness of the edges
        real(wp) :: edgetan, initial, v_norm, xnormal, ynormal !< The manifoldness of the edges

        ! Total number of edges in 2D STL
        edge_count = 3*model%ntrs

        ! Initialize edge_occurrence array to zero
        edge_occurrence = 0
        edge_index = 0

        ! Collect all edges of all triangles and store them
        do i = 1, model%ntrs
            ! First edge (v1, v2)
            edge(1, 1:2) = model%trs(i)%v(1, 1:2)
            edge(2, 1:2) = model%trs(i)%v(2, 1:2)
            call f_register_edge(temp_boundary_v, edge, edge_index, edge_count)

            ! Second edge (v2, v3)
            edge(1, 1:2) = model%trs(i)%v(2, 1:2)
            edge(2, 1:2) = model%trs(i)%v(3, 1:2)
            call f_register_edge(temp_boundary_v, edge, edge_index, edge_count)

            ! Third edge (v3, v1)
            edge(1, 1:2) = model%trs(i)%v(3, 1:2)
            edge(2, 1:2) = model%trs(i)%v(1, 1:2)
            call f_register_edge(temp_boundary_v, edge, edge_index, edge_count)
        end do

        ! Check all edges and count repeated edges
        $:GPU_PARALLEL_LOOP(private='[i,j]', copy='[temp_boundary_v,edge_occurrence]', copyin='[threshold_edge_zero]', collapse=2)
        do i = 1, edge_count
            do j = 1, edge_count
                if (i /= j) then
                    if (((abs(temp_boundary_v(i, 1, 1) - temp_boundary_v(j, 1, 1)) < threshold_edge_zero) .and. &
                         (abs(temp_boundary_v(i, 1, 2) - temp_boundary_v(j, 1, 2)) < threshold_edge_zero) .and. &
                         (abs(temp_boundary_v(i, 2, 1) - temp_boundary_v(j, 2, 1)) < threshold_edge_zero) .and. &
                         (abs(temp_boundary_v(i, 2, 2) - temp_boundary_v(j, 2, 2)) < threshold_edge_zero)) .or. &
                        ((abs(temp_boundary_v(i, 1, 1) - temp_boundary_v(j, 2, 1)) < threshold_edge_zero) .and. &
                         (abs(temp_boundary_v(i, 1, 2) - temp_boundary_v(j, 2, 2)) < threshold_edge_zero) .and. &
                         (abs(temp_boundary_v(i, 2, 1) - temp_boundary_v(j, 1, 1)) < threshold_edge_zero) .and. &
                         (abs(temp_boundary_v(i, 2, 2) - temp_boundary_v(j, 1, 2)) < threshold_edge_zero))) then

                        $:GPU_ATOMIC(atomic='update')
                        edge_occurrence(i) = edge_occurrence(i) + 1
                    end if
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        ! Count the number of boundary vertices/edges
        boundary_vertex_count = 0
        boundary_edge_count = 0

        do i = 1, edge_count
            if (edge_occurrence(i) == 0) then
                boundary_vertex_count = boundary_vertex_count + 2
                boundary_edge_count = boundary_edge_count + 1
            end if
        end do

        ! Allocate the boundary_v array based on the number of boundary edges
        allocate (boundary_v(boundary_edge_count, 1:3, 1:2))

        ! Store boundary vertices
        store_index = 0
        do i = 1, edge_count
            if (edge_occurrence(i) == 0) then
                store_index = store_index + 1
                boundary_v(store_index, 1, 1:2) = temp_boundary_v(i, 1, 1:2)
                boundary_v(store_index, 2, 1:2) = temp_boundary_v(i, 2, 1:2)
            end if
        end do

        ! Find/store the normal vector of the boundary edges
        do i = 1, boundary_edge_count
            boundary_edge(1) = boundary_v(i, 2, 1) - boundary_v(i, 1, 1)
            boundary_edge(2) = boundary_v(i, 2, 2) - boundary_v(i, 1, 2)
            edgetan = boundary_edge(1)/sign(max(sgm_eps, abs(boundary_edge(2))), boundary_edge(2))

            if (abs(boundary_edge(2)) < threshold_vector_zero) then
                if (edgetan > 0._wp) then
                    ynormal = -1
                    xnormal = 0._wp
                else
                    ynormal = 1
                    xnormal = 0._wp
                end if
            else
                initial = boundary_edge(2)
                ynormal = -edgetan*initial
                xnormal = initial
            end if

            v_norm = sqrt(xnormal**2 + ynormal**2)
            boundary_v(i, 3, 1) = xnormal/v_norm
            boundary_v(i, 3, 2) = ynormal/v_norm
        end do

    end subroutine f_check_boundary

    !> This procedure appends the edge end vertices to a temporary buffer.
    subroutine f_register_edge(temp_boundary_v, edge, edge_index, edge_count)

        integer, intent(inout) :: edge_index !< Edge index iterator
        integer, intent(inout) :: edge_count !< Total number of edges
        real(wp), intent(in), dimension(1:2, 1:2) :: edge !< Edges end points to be registered
        real(wp), dimension(1:edge_count, 1:2, 1:2), intent(inout) :: temp_boundary_v !< Temporary edge end vertex buffer

        ! Increment edge index and store the edge
        edge_index = edge_index + 1
        temp_boundary_v(edge_index, 1, 1:2) = edge(1, 1:2)
        temp_boundary_v(edge_index, 2, 1:2) = edge(2, 1:2)

    end subroutine f_register_edge

    !> This procedure determines the levelset distance and normals of 3D models
    !! by computing the exact closest point via projection onto triangle surfaces.
    !! @param ntrs                  Number of triangles for this patch
    !! @param trs_v                 Flat GPU array of triangle vertices for all patches
    !! @param trs_n                 Flat GPU array of triangle normals for all patches
    !! @param pid                   Patch index into the arrays
    !! @param point                 The cell center of the current levelset cell
    !! @param normals               Output levelset normals
    !! @param distance              Output levelset distance
    subroutine f_distance_normals_3D(ntrs, trs_v, trs_n, pid, point, normals, distance)
        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in) :: ntrs
        real(wp), dimension(:, :, :, :), intent(in) :: trs_v
        real(wp), dimension(:, :, :), intent(in) :: trs_n
        integer, intent(in) :: pid
        real(wp), dimension(1:3), intent(in) :: point
        real(wp), dimension(1:3), intent(out) :: normals
        real(wp), intent(out) :: distance

        integer :: i, j
        real(wp) :: dist_min, dist_proj, dist_v, dist_e, t
        real(wp) :: v1(1:3), v2(1:3), v3(1:3)
        real(wp) :: e0(1:3), e1(1:3), pv(1:3)
        real(wp) :: n(1:3), proj(1:3)
        real(wp) :: d, ndot, denom
        real(wp) :: u, v_bary, w
        real(wp) :: l00, l01, l11, l20, l21
        real(wp) :: edge(1:3), pe(1:3)
        real(wp) :: verts(1:3, 1:3)

        dist_min = 1.e12_wp
        normals = 0._wp

        do i = 1, ntrs
            ! Triangle vertices
            v1(:) = trs_v(1, :, i, pid)
            v2(:) = trs_v(2, :, i, pid)
            v3(:) = trs_v(3, :, i, pid)

            ! Triangle normal
            n(:) = trs_n(:, i, pid)

            ! Project point onto triangle plane
            pv(:) = point(:) - v1(:)
            d = dot_product(pv, n)
            if (abs(d) >= dist_min) cycle ! minimum distance is not small enough, no need to check validity
            proj(:) = point(:) - d*n(:)

            ! Check if projection is inside triangle using barycentric coordinates
            e0(:) = v2(:) - v1(:)
            e1(:) = v3(:) - v1(:)
            pv(:) = proj(:) - v1(:)

            l00 = dot_product(e0, e0)
            l01 = dot_product(e0, e1)
            l11 = dot_product(e1, e1)
            l20 = dot_product(pv, e0)
            l21 = dot_product(pv, e1)

            denom = l00*l11 - l01*l01

            if (abs(denom) > 0._wp) then
                v_bary = (l11*l20 - l01*l21)/denom
                w = (l00*l21 - l01*l20)/denom
                u = 1._wp - v_bary - w
            else
                u = -1._wp
                v_bary = -1._wp
                w = -1._wp
            end if

            ! If projection is inside triangle
            if (u >= 0._wp .and. v_bary >= 0._wp .and. w >= 0._wp) then
                dist_proj = sqrt((point(1) - proj(1))**2 + &
                                 (point(2) - proj(2))**2 + &
                                 (point(3) - proj(3))**2)

                if (dist_proj < dist_min) then
                    dist_min = dist_proj
                    normals(:) = n(:)
                end if
            else
                ! Projection outside triangle: check edges and vertices
                verts(:, 1) = v1(:)
                verts(:, 2) = v2(:)
                verts(:, 3) = v3(:)

                ! Check three edges
                do j = 1, 3
                    edge(:) = verts(:, mod(j, 3) + 1) - verts(:, j)
                    pe(:) = point(:) - verts(:, j)

                    t = dot_product(pe, edge)/max(dot_product(edge, edge), 1.e-30_wp)

                    if (t >= 0._wp .and. t <= 1._wp) then
                        proj(:) = verts(:, j) + t*edge(:)
                        dist_e = sqrt((point(1) - proj(1))**2 + &
                                      (point(2) - proj(2))**2 + &
                                      (point(3) - proj(3))**2)

                        if (dist_e < dist_min) then
                            dist_min = dist_e
                            normals(:) = n(:)
                        end if
                    end if
                end do

                ! Check three vertices
                do j = 1, 3
                    dist_v = sqrt((point(1) - verts(1, j))**2 + &
                                  (point(2) - verts(2, j))**2 + &
                                  (point(3) - verts(3, j))**2)

                    if (dist_v < dist_min) then
                        dist_min = dist_v
                        normals(:) = n(:)
                    end if
                end do
            end if
        end do

        distance = dist_min

    end subroutine f_distance_normals_3D

    !> This procedure determines the levelset distance and normals of 2D models
    !! by computing the exact closest point via projection onto boundary edges.
    !! @param boundary_v            Flat GPU array of boundary vertices/normals for all patches
    !! @param pid                   Patch index into the boundary_v array
    !! @param boundary_edge_count   Total number of boundary edges for this patch
    !! @param point                 The cell center of the current levelset cell
    !! @param normals               Output levelset normals
    !! @param distance              Output levelset distance
    subroutine f_distance_normals_2D(boundary_v, pid, boundary_edge_count, point, normals, distance)
        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(:, :, :, :), intent(in) :: boundary_v
        integer, intent(in) :: pid
        integer, intent(in) :: boundary_edge_count
        real(wp), dimension(1:3), intent(in) :: point
        real(wp), dimension(1:3), intent(out) :: normals
        real(wp), intent(out) :: distance

        integer :: i
        real(wp) :: dist_min, dist, t
        real(wp) :: v1(1:2), v2(1:2), edge(1:2), pv(1:2)
        real(wp) :: edge_len_sq, proj(1:2)

        dist_min = 1.e12_wp
        normals = 0._wp

        do i = 1, boundary_edge_count
            ! Edge endpoints
            v1(1) = boundary_v(i, 1, 1, pid)
            v1(2) = boundary_v(i, 1, 2, pid)
            v2(1) = boundary_v(i, 2, 1, pid)
            v2(2) = boundary_v(i, 2, 2, pid)

            ! Edge vector and point-to-v1 vector
            edge = v2 - v1
            pv(1) = point(1) - v1(1)
            pv(2) = point(2) - v1(2)
            edge_len_sq = dot_product(edge, edge)

            ! Parameter of projection onto the edge line
            if (edge_len_sq > 0._wp) then
                t = dot_product(pv, edge)/edge_len_sq
            else
                t = 0._wp
            end if

            ! Check if projection falls within the segment
            if (t >= 0._wp .and. t <= 1._wp) then
                proj = v1 + t*edge
                dist = sqrt((point(1) - proj(1))**2 + (point(2) - proj(2))**2)
            else if (t < 0._wp) then ! negative t means that v1 is the closest point on the edge
                dist = sqrt((point(1) - v1(1))**2 + (point(2) - v1(2))**2)
            else ! t > 1 means that v2 is the closest point on the line edge
                dist = sqrt((point(1) - v2(1))**2 + (point(2) - v2(2))**2)
            end if

            if (dist < dist_min) then
                dist_min = dist
                normals(1) = boundary_v(i, 3, 1, pid)
                normals(2) = boundary_v(i, 3, 2, pid)
            end if
        end do

        distance = dist_min

    end subroutine f_distance_normals_2D

    subroutine s_pack_model_for_gpu(ma)
        type(t_model_array), intent(inout) :: ma
        integer :: i

        ma%ntrs = ma%model%ntrs
        allocate (ma%trs_v(1:3, 1:3, 1:ma%ntrs))
        allocate (ma%trs_n(1:3, 1:ma%ntrs))

        do i = 1, ma%ntrs
            ma%trs_v(:, :, i) = ma%model%trs(i)%v(:, :)
            ma%trs_n(:, i) = ma%model%trs(i)%n(:)
        end do
    end subroutine

end module m_model
