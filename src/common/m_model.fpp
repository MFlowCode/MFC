!>
!! @file   m_model.fpp
!! @author Henry Le Berre <hberre3@gatech.edu>
!! @brief  Contains module m_model

#:include 'macros.fpp'

module m_model

    use m_helper
    use m_mpi_proxy
    use m_derived_types

    use iso_c_binding, only: c_char, c_int32_t, c_int16_t, c_float

    implicit none

    private

    public :: f_model_read, s_model_write, s_model_free, f_model_is_inside

    ! Subroutines for STL immersed boundaries
    public :: f_check_boundary, f_register_edge, f_check_interpolation_2D, &
              f_check_interpolation_3D, f_interpolate_2D, f_interpolate_3D, &
              f_interpolated_distance, f_normals, f_distance, f_distance_normals_3D, f_tri_area

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

        real(kind=c_float) :: normal(3), v(3, 3)
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
            read (line(13:), *) model%trs(i)%n

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

    !> This procedure, recursively, finds whether a point is inside an octree.
    !! @param model    Model to search in.
    !! @param point    Point to test.
    !! @param spacing  Space around the point to search in (grid spacing).
    !! @param spc      Number of samples per cell.
    !! @return True if the point is inside the octree, false otherwise.
    impure function f_model_is_inside(model, point, spacing, spc) result(fraction)

        type(t_model), intent(in) :: model
        real(wp), dimension(1:3), intent(in) :: point
        real(wp), dimension(1:3), intent(in) :: spacing
        integer, intent(in) :: spc

        real(wp) :: fraction

        type(t_ray) :: ray
        integer :: i, j, nInOrOut, nHits

        real(wp), dimension(1:spc, 1:3) :: ray_origins, ray_dirs

        do i = 1, spc
            call random_number(ray_origins(i, :))
            ray_origins(i, :) = point + (ray_origins(i, :) - 0.5_wp)*spacing(:)

            call random_number(ray_dirs(i, :))
            ray_dirs(i, :) = ray_dirs(i, :) - 0.5_wp
            ray_dirs(i, :) = ray_dirs(i, :)/sqrt(sum(ray_dirs(i, :)*ray_dirs(i, :)))
        end do

        nInOrOut = 0
        do i = 1, spc
            ray%o = ray_origins(i, :)
            ray%d = ray_dirs(i, :)

            nHits = 0
            do j = 1, model%ntrs
                if (f_intersects_triangle(ray, model%trs(j))) then
                    nHits = nHits + 1
                end if
            end do

            nInOrOut = nInOrOut + mod(nHits, 2)
        end do

        fraction = real(nInOrOut)/real(spc)

    end function f_model_is_inside

    ! From https://www.scratchapixel.com/lessons/3e-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
    !> This procedure checks if a ray intersects a triangle.
    !! @param ray      Ray.
    !! @param triangle Triangle.
    !! @return         True if the ray intersects the triangle, false otherwise.
    elemental function f_intersects_triangle(ray, triangle) result(intersects)

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
    !! @param boundary_v                 Output boundary vertices/normals.
    !! @param boundary_vertex_count      Output total boundary vertex count
    !! @param boundary_edge_count        Output total boundary edge counts
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

                        edge_occurrence(i) = edge_occurrence(i) + 1
                    end if
                end if
            end do
        end do

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
            edgetan = boundary_edge(1)/boundary_edge(2)

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
    !! @param temp_boundary_v      Temporary edge end vertex buffer
    !! @param edge                 Edges end points to be registered
    !! @param edge_index           Edge index iterator
    !! @param edge_count           Total number of edges
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

    !> This procedure check if interpolates is needed for 2D models.
    !! @param boundary_v                Temporary edge end vertex buffer
    !! @param boundary_edge_count       Output total number of boundary edges
    !! @param spacing                   Dimensions of the current levelset cell
    !! @param interpolate               Logical output
    subroutine f_check_interpolation_2D(boundary_v, boundary_edge_count, spacing, interpolate)

        logical, intent(inout) :: interpolate !< Logical indicator of interpolation
        integer, intent(in) :: boundary_edge_count !< Number of boundary edges
        real(wp), intent(in), dimension(1:boundary_edge_count, 1:3, 1:2) :: boundary_v
        real(wp), dimension(1:3), intent(in) :: spacing

        real(wp) :: l1, cell_width !< Length of each boundary edge and cell width
        integer :: j !< Boundary edge index iterator

        cell_width = minval(spacing(1:2))
        interpolate = .false.

        do j = 1, boundary_edge_count

            l1 = sqrt((boundary_v(j, 2, 1) - boundary_v(j, 1, 1))**2 + &
                      (boundary_v(j, 2, 2) - boundary_v(j, 1, 2))**2)

            if ((l1 > cell_width)) then
                interpolate = .true.
            else
                interpolate = .false.
            end if
        end do

    end subroutine f_check_interpolation_2D

    !> This procedure check if interpolates is needed for 3D models.
    !! @param model              Model to search in.
    !! @param spacing            Dimensions of the current levelset cell
    !! @param interpolate        Logical output
    subroutine f_check_interpolation_3D(model, spacing, interpolate)

        logical, intent(inout) :: interpolate
        type(t_model), intent(in) :: model
        real(wp), dimension(1:3), intent(in) :: spacing
        real(wp), dimension(1:3) :: edge_l
        real(wp) :: cell_width
        real(wp), dimension(1:3, 1:3) :: tri_v
        integer :: i, j !< Loop iterator

        cell_width = minval(spacing)
        interpolate = .false.

        do i = 1, model%ntrs
            do j = 1, 3
                tri_v(1, j) = model%trs(i)%v(1, j)
                tri_v(2, j) = model%trs(i)%v(2, j)
                tri_v(3, j) = model%trs(i)%v(3, j)
            end do

            edge_l(1) = sqrt((tri_v(1, 2) - tri_v(1, 1))**2 + &
                             (tri_v(2, 2) - tri_v(2, 1))**2 + &
                             (tri_v(3, 2) - tri_v(3, 1))**2)
            edge_l(2) = sqrt((tri_v(1, 3) - tri_v(1, 2))**2 + &
                             (tri_v(2, 3) - tri_v(2, 2))**2 + &
                             (tri_v(3, 3) - tri_v(3, 2))**2)
            edge_l(3) = sqrt((tri_v(1, 1) - tri_v(1, 3))**2 + &
                             (tri_v(2, 1) - tri_v(2, 3))**2 + &
                             (tri_v(3, 1) - tri_v(3, 3))**2)

            if ((edge_l(1) > cell_width) .or. &
                (edge_l(2) > cell_width) .or. &
                (edge_l(3) > cell_width)) then
                interpolate = .true.
            else
                interpolate = .false.
            end if
        end do

    end subroutine f_check_interpolation_3D

    !> This procedure interpolates 2D models.
    !! @param boundary_v                   Group of all the boundary vertices of the 2D model without interpolation
    !! @param boundary_edge_count          Output total number of boundary edges
    !! @param spacing                      Dimensions of the current levelset cell
    !! @param interpolated_boundary_v      Output all the boundary vertices of the interpolated 2D model
    !! @param total_vertices               Total number of vertices after interpolation
    subroutine f_interpolate_2D(boundary_v, boundary_edge_count, spacing, interpolated_boundary_v, total_vertices)

        real(wp), intent(in), dimension(:, :, :) :: boundary_v
        real(wp), dimension(1:3), intent(in) :: spacing
        real(wp), allocatable, intent(inout), dimension(:, :) :: interpolated_boundary_v

        integer, intent(inout) :: total_vertices, boundary_edge_count
        integer :: num_segments
        integer :: i, j

        real(wp) :: edge_length, cell_width
        real(wp), dimension(1:2) :: edge_x, edge_y, edge_del

        ! Get the number of boundary edges
        cell_width = minval(spacing(1:2))
        num_segments = 0

        ! First pass: Calculate the total number of vertices including interpolated ones
        total_vertices = 1
        do i = 1, boundary_edge_count
            ! Get the coordinates of the two ends of the current edge
            edge_x(1) = boundary_v(i, 1, 1)
            edge_y(1) = boundary_v(i, 1, 2)
            edge_x(2) = boundary_v(i, 2, 1)
            edge_y(2) = boundary_v(i, 2, 2)

            ! Compute the length of the edge
            edge_length = sqrt((edge_x(2) - edge_x(1))**2 + &
                               (edge_y(2) - edge_y(1))**2)

            ! Determine the number of segments
            if (edge_length > cell_width) then
                num_segments = Ifactor_2D*ceiling(edge_length/cell_width)
            else
                num_segments = 1
            end if

            ! Each edge contributes num_segments vertices
            total_vertices = total_vertices + num_segments
        end do

        ! Allocate memory for the new boundary vertices array
        allocate (interpolated_boundary_v(1:total_vertices, 1:3))

        ! Fill the new boundary vertices array with original and interpolated vertices
        total_vertices = 1
        do i = 1, boundary_edge_count
            ! Get the coordinates of the two ends of the current edge
            edge_x(1) = boundary_v(i, 1, 1)
            edge_y(1) = boundary_v(i, 1, 2)
            edge_x(2) = boundary_v(i, 2, 1)
            edge_y(2) = boundary_v(i, 2, 2)

            ! Compute the length of the edge
            edge_length = sqrt((edge_x(2) - edge_x(1))**2 + &
                               (edge_y(2) - edge_y(1))**2)

            ! Determine the number of segments and interpolation step
            if (edge_length > cell_width) then
                num_segments = Ifactor_2D*ceiling(edge_length/cell_width)
                edge_del(1) = (edge_x(2) - edge_x(1))/num_segments
                edge_del(2) = (edge_y(2) - edge_y(1))/num_segments
            else
                num_segments = 1
                edge_del(1) = 0._wp
                edge_del(2) = 0._wp
            end if

            interpolated_boundary_v(1, 1) = edge_x(1)
            interpolated_boundary_v(1, 2) = edge_y(1)
            interpolated_boundary_v(1, 3) = 0._wp

            ! Add original and interpolated vertices to the output array
            do j = 1, num_segments - 1
                total_vertices = total_vertices + 1
                interpolated_boundary_v(total_vertices, 1) = edge_x(1) + j*edge_del(1)
                interpolated_boundary_v(total_vertices, 2) = edge_y(1) + j*edge_del(2)
            end do

            ! Add the last vertex of the edge
            if (num_segments > 0) then
                total_vertices = total_vertices + 1
                interpolated_boundary_v(total_vertices, 1) = edge_x(2)
                interpolated_boundary_v(total_vertices, 2) = edge_y(2)
            end if
        end do

    end subroutine f_interpolate_2D

    !> This procedure interpolates 3D models.
    !! @param model                        Model to search in.
    !! @param spacing                      Dimensions of the current levelset cell
    !! @param interpolated_boundary_v      Output all the boundary vertices of the interpolated 3D model
    !! @param total_vertices               Total number of vertices after interpolation
    impure subroutine f_interpolate_3D(model, spacing, interpolated_boundary_v, total_vertices)
        real(wp), dimension(1:3), intent(in) :: spacing
        type(t_model), intent(in) :: model
        real(wp), allocatable, intent(inout), dimension(:, :) :: interpolated_boundary_v
        integer, intent(out) :: total_vertices

        integer :: i, j, k, num_triangles, num_segments, num_inner_vertices
        real(wp), dimension(1:3, 1:3) :: tri
        real(wp), dimension(1:3) :: edge_del, cell_area
        real(wp), dimension(1:3) :: bary_coord !< Barycentric coordinates
        real(wp) :: edge_length, cell_width, cell_area_min, tri_area

        ! Number of triangles in the model
        num_triangles = model%ntrs
        cell_width = minval(spacing)

        ! Find the minimum surface area
        cell_area(1) = spacing(1)*spacing(2)
        cell_area(2) = spacing(1)*spacing(3)
        cell_area(3) = spacing(2)*spacing(3)
        cell_area_min = minval(cell_area)
        num_inner_vertices = 0

        ! Calculate the total number of vertices including interpolated ones
        total_vertices = 0
        do i = 1, num_triangles
            do j = 1, 3
                ! Get the coordinates of the two vertices of the current edge
                tri(1, 1) = model%trs(i)%v(j, 1)
                tri(1, 2) = model%trs(i)%v(j, 2)
                tri(1, 3) = model%trs(i)%v(j, 3)
                ! Next vertex in the triangle (cyclic)
                tri(2, 1) = model%trs(i)%v(mod(j, 3) + 1, 1)
                tri(2, 2) = model%trs(i)%v(mod(j, 3) + 1, 2)
                tri(2, 3) = model%trs(i)%v(mod(j, 3) + 1, 3)

                ! Compute the length of the edge
                edge_length = sqrt((tri(2, 1) - tri(1, 1))**2 + &
                                   (tri(2, 2) - tri(1, 2))**2 + &
                                   (tri(2, 3) - tri(1, 3))**2)

                ! Determine the number of segments
                if (edge_length > cell_width) then
                    num_segments = Ifactor_3D*ceiling(edge_length/cell_width)
                else
                    num_segments = 1
                end if

                ! Each edge contributes num_segments vertices
                total_vertices = total_vertices + num_segments + 1
            end do

            ! Add vertices inside the triangle
            do k = 1, 3
                tri(k, 1) = model%trs(i)%v(k, 1)
                tri(k, 2) = model%trs(i)%v(k, 2)
                tri(k, 3) = model%trs(i)%v(k, 3)
            end do
            call f_tri_area(tri, tri_area)

            if (tri_area > threshold_bary*cell_area_min) then
                num_inner_vertices = Ifactor_bary_3D*ceiling(tri_area/cell_area_min)
                total_vertices = total_vertices + num_inner_vertices
            end if
        end do

        ! Allocate memory for the new boundary vertices array
        allocate (interpolated_boundary_v(1:total_vertices, 1:3))

        ! Fill the new boundary vertices array with original and interpolated vertices
        total_vertices = 0
        do i = 1, num_triangles
            ! Loop through the 3 edges of each triangle
            do j = 1, 3
                ! Get the coordinates of the two vertices of the current edge
                tri(1, 1) = model%trs(i)%v(j, 1)
                tri(1, 2) = model%trs(i)%v(j, 2)
                tri(1, 3) = model%trs(i)%v(j, 3)
                ! Next vertex in the triangle (cyclic)
                tri(2, 1) = model%trs(i)%v(mod(j, 3) + 1, 1)
                tri(2, 2) = model%trs(i)%v(mod(j, 3) + 1, 2)
                tri(2, 3) = model%trs(i)%v(mod(j, 3) + 1, 3)

                ! Compute the length of the edge
                edge_length = sqrt((tri(2, 1) - tri(1, 1))**2 + &
                                   (tri(2, 2) - tri(1, 2))**2 + &
                                   (tri(2, 3) - tri(1, 3))**2)

                ! Determine the number of segments and interpolation step
                if (edge_length > cell_width) then
                    num_segments = Ifactor_3D*ceiling(edge_length/cell_width)
                    edge_del(1) = (tri(2, 1) - tri(1, 1))/num_segments
                    edge_del(2) = (tri(2, 2) - tri(1, 2))/num_segments
                    edge_del(3) = (tri(2, 3) - tri(1, 3))/num_segments
                else
                    num_segments = 1
                    edge_del = 0._wp
                end if

                ! Add original and interpolated vertices to the output array
                do k = 0, num_segments - 1
                    total_vertices = total_vertices + 1
                    interpolated_boundary_v(total_vertices, 1) = tri(1, 1) + k*edge_del(1)
                    interpolated_boundary_v(total_vertices, 2) = tri(1, 2) + k*edge_del(2)
                    interpolated_boundary_v(total_vertices, 3) = tri(1, 3) + k*edge_del(3)
                end do

                ! Add the last vertex of the edge
                total_vertices = total_vertices + 1
                interpolated_boundary_v(total_vertices, 1) = tri(2, 1)
                interpolated_boundary_v(total_vertices, 2) = tri(2, 2)
                interpolated_boundary_v(total_vertices, 3) = tri(2, 3)
            end do

            ! Interpolate verties that are not on edges
            do k = 1, 3
                tri(k, 1) = model%trs(i)%v(k, 1)
                tri(k, 2) = model%trs(i)%v(k, 2)
                tri(k, 3) = model%trs(i)%v(k, 3)
            end do
            call f_tri_area(tri, tri_area)

            if (tri_area > threshold_bary*cell_area_min) then
                num_inner_vertices = Ifactor_bary_3D*ceiling(tri_area/cell_area_min)
                !Use barycentric coordinates for randomly distributed points
                do k = 1, num_inner_vertices
                    call random_number(bary_coord(1))
                    call random_number(bary_coord(2))

                    if ((bary_coord(1) + bary_coord(2)) >= 1._wp) then
                        bary_coord(1) = 1._wp - bary_coord(1)
                        bary_coord(2) = 1._wp - bary_coord(2)
                    end if
                    bary_coord(3) = 1._wp - bary_coord(1) - bary_coord(2)

                    total_vertices = total_vertices + 1
                    interpolated_boundary_v(total_vertices, 1) = dot_product(bary_coord, tri(1:3, 1))
                    interpolated_boundary_v(total_vertices, 2) = dot_product(bary_coord, tri(1:3, 2))
                    interpolated_boundary_v(total_vertices, 3) = dot_product(bary_coord, tri(1:3, 3))
                end do
            end if
        end do

    end subroutine f_interpolate_3D

    !> This procedure determines the levelset distance and normals of the 3D models without interpolation.
    !! @param model        Model to search in.
    !! @param point        The cell centers of the current level cell
    !! @param normals      The output levelset normals
    !! @param distance     The output levelset distance
    subroutine f_distance_normals_3D(model, point, normals, distance)

        type(t_model), intent(IN) :: model
        real(wp), dimension(1:3), intent(in) :: point
        real(wp), dimension(1:3), intent(out) :: normals
        real(wp), intent(out) :: distance

        real(wp), dimension(1:3, 1:3) :: tri
        real(wp) :: dist_min, dist_t_min
        real(wp) :: dist_min_normal, dist_buffer_normal
        real(wp), dimension(1:3) :: midp !< Centers of the triangle facets
        real(wp), dimension(1:3) :: dist_buffer !< Distance between the cell center and the vertices
        integer :: i, j, tri_idx !< Iterator

        dist_min = 1.e12_wp
        dist_min_normal = 1.e12_wp
        distance = 0._wp

        tri_idx = 0
        do i = 1, model%ntrs
            do j = 1, 3
                tri(j, 1) = model%trs(i)%v(j, 1)
                tri(j, 2) = model%trs(i)%v(j, 2)
                tri(j, 3) = model%trs(i)%v(j, 3)
                dist_buffer(j) = sqrt((point(1) - tri(j, 1))**2 + &
                                      (point(2) - tri(j, 2))**2 + &
                                      (point(3) - tri(j, 3))**2)
            end do

            ! Get the surface center of each triangle facet
            do j = 1, 3
                midp(j) = (tri(1, j) + tri(2, j) + tri(3, j))/3
            end do

            dist_t_min = minval(dist_buffer(1:3))
            dist_buffer_normal = sqrt((point(1) - midp(1))**2 + &
                                      (point(2) - midp(2))**2 + &
                                      (point(3) - midp(3))**2)

            if (dist_t_min < dist_min) then
                dist_min = dist_t_min
            end if

            if (dist_buffer_normal < dist_min_normal) then
                dist_min_normal = dist_buffer_normal
                tri_idx = i
            end if
        end do

        normals(1) = model%trs(tri_idx)%n(1)
        normals(2) = model%trs(tri_idx)%n(2)
        normals(3) = model%trs(tri_idx)%n(3)
        distance = dist_min

    end subroutine f_distance_normals_3D

    !> This procedure determines the levelset distance of 2D models without interpolation.
    !! @param boundary_v                   Group of all the boundary vertices of the 2D model without interpolation
    !! @param boundary_vertex_count        Output the total number of boundary vertices
    !! @param boundary_edge_count          Output the total number of boundary edges
    !! @param point                        The cell centers of the current levelset cell
    !! @param spacing                      Dimensions of the current levelset cell
    !! @return                             Distance which the levelset distance without interpolation
    function f_distance(boundary_v, boundary_edge_count, point) result(distance)

        integer, intent(in) :: boundary_edge_count
        real(wp), intent(in), dimension(1:boundary_edge_count, 1:3, 1:2) :: boundary_v
        real(wp), dimension(1:3), intent(in) :: point

        integer :: i
        real(wp) :: dist_buffer1, dist_buffer2
        real(wp), dimension(1:boundary_edge_count) :: dist_buffer
        real(wp) :: distance

        distance = 0._wp
        do i = 1, boundary_edge_count
            dist_buffer1 = sqrt((point(1) - boundary_v(i, 1, 1))**2 + &
                                & (point(2) - boundary_v(i, 1, 2))**2)

            dist_buffer2 = sqrt((point(1) - boundary_v(i, 2, 1))**2 + &
                                & (point(2) - boundary_v(i, 2, 2))**2)

            dist_buffer(i) = minval((/dist_buffer1, dist_buffer2/))
        end do

        distance = minval(dist_buffer)

    end function f_distance

    !> This procedure determines the levelset normals of 2D models without interpolation.
    !! @param boundary_v                   Group of all the boundary vertices of the 2D model without interpolation
    !! @param boundary_edge_count          Output the total number of boundary edges
    !! @param point                        The cell centers of the current levelset cell
    !! @param normals                      Output levelset normals without interpolation
    subroutine f_normals(boundary_v, boundary_edge_count, point, normals)

        integer, intent(in) :: boundary_edge_count
        real(wp), intent(in), dimension(1:boundary_edge_count, 1:3, 1:2) :: boundary_v
        real(wp), dimension(1:3), intent(in) :: point
        real(wp), dimension(1:3), intent(out) :: normals

        integer :: i, idx_buffer
        real(wp) :: dist_min, dist_buffer
        real(wp) :: midp(1:3)

        dist_buffer = 0._wp
        dist_min = initial_distance_buffer
        idx_buffer = 0

        do i = 1, boundary_edge_count
            midp(1) = (boundary_v(i, 2, 1) + boundary_v(i, 1, 1))/2
            midp(2) = (boundary_v(i, 2, 2) + boundary_v(i, 1, 2))/2
            midp(3) = 0._wp

            dist_buffer = sqrt((point(1) - midp(1))**2 + &
                                & (point(2) - midp(2))**2)

            if (dist_buffer < dist_min) then
                dist_min = dist_buffer
                idx_buffer = i
            end if
        end do

        normals(1) = boundary_v(idx_buffer, 3, 1)
        normals(2) = boundary_v(idx_buffer, 3, 2)
        normals(3) = 0._wp

    end subroutine f_normals

    !> This procedure calculates the barycentric facet area
    subroutine f_tri_area(tri, tri_area)

        real(wp), dimension(1:3, 1:3), intent(in) :: tri
        real(wp), intent(out) :: tri_area
        real(wp), dimension(1:3) :: AB, AC, cross
        integer :: i !< Loop iterator

        do i = 1, 3
            AB(i) = tri(2, i) - tri(1, i)
            AC(i) = tri(3, i) - tri(1, i)
        end do

        cross(1) = AB(2)*AC(3) - AB(3)*AC(2)
        cross(2) = AB(3)*AC(1) - AB(1)*AC(3)
        cross(3) = AB(1)*AC(2) - AB(2)*AC(1)
        tri_area = 0.5_wp*sqrt(cross(1)**2 + cross(2)**2 + cross(3)**2)

    end subroutine f_tri_area

    !> This procedure determines the levelset of interpolated 2D models.
    !! @param interpolated_boundary_v      Group of all the boundary vertices of the interpolated 2D model
    !! @param total_vertices               Total number of vertices after interpolation
    !! @param point                        The cell centers of the current levelset cell
    !! @return                             Distance which the levelset distance without interpolation
    function f_interpolated_distance(interpolated_boundary_v, total_vertices, point) result(distance)

        integer, intent(in) :: total_vertices
        real(wp), intent(in), dimension(1:total_vertices, 1:3) :: interpolated_boundary_v
        real(wp), dimension(1:3), intent(in) :: point

        integer :: i !< Loop iterator
        real(wp) :: dist_buffer, min_dist
        real(wp) :: distance

        distance = initial_distance_buffer
        dist_buffer = initial_distance_buffer
        min_dist = initial_distance_buffer

        do i = 1, total_vertices
            dist_buffer = sqrt((point(1) - interpolated_boundary_v(i, 1))**2 + &
                               (point(2) - interpolated_boundary_v(i, 2))**2 + &
                               (point(3) - interpolated_boundary_v(i, 3))**2)

            if (min_dist > dist_buffer) then
                min_dist = dist_buffer
            end if
        end do

        distance = min_dist

    end function f_interpolated_distance

end module m_model
