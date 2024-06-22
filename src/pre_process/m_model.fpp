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

contains

    !> This procedure reads a binary STL file.
    !! @param filepath Path to the STL file.
    !! @param model The binary of the STL file.
    subroutine s_read_stl_binary(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(out) :: model

        integer :: i, j, iunit, iostat

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
    subroutine s_read_stl_ascii(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(out) :: model

        integer :: i, j, iunit, iostat

        character(80) :: line

        open (newunit=iunit, file=filepath, action='READ', &
              form='FORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open ASCII STL file ", filepath

            call s_mpi_abort()
        end if

        model%ntrs = 0
        do
            if (.not. f_read_line(iunit, line)) exit

            if (line(1:6) == "facet ") then
                model%ntrs = model%ntrs + 1
            end if
        end do

        allocate (model%trs(model%ntrs))

        rewind (iunit)

        i = 1
        do
            if (.not. f_read_line(iunit, line)) exit

            if (line(1:5) == "solid") cycle
            if (line(1:8) == "endsolid") exit

            if (line(1:12) /= "facet normal") then
                print *, "Error: expected facet normal in STL file ", filepath

                call s_mpi_abort()
            end if

            call s_skip_ignored_lines(iunit)
            read (line(13:), *) model%trs(i)%n

            call s_skip_ignored_lines(iunit)
            read (iunit, '(A)') line

            do j = 1, 3
                if (.not. f_read_line(iunit, line)) exit

                if (line(1:6) /= "vertex") then
                    print *, "Error: expected vertex in STL file ", filepath

                    call s_mpi_abort()
                end if

                call s_skip_ignored_lines(iunit)
                read (line(7:), *) model%trs(i)%v(j, :)
            end do

            if (.not. f_read_line(iunit, line)) exit
            if (.not. f_read_line(iunit, line)) exit

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
    subroutine s_read_stl(filepath, model)

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

    end subroutine

    !> This procedure reads an OBJ file.
    !! @param filepath Path to the odj file.
    !! @param model The obj file.
    subroutine s_read_obj(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(out) :: model

        integer :: i, j, k, l, iunit, iostat, nVertices

        t_vec3, allocatable :: vertices(:, :)

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

    end subroutine

    !> This procedure reads a mesh from a file.
    !! @param filepath Path to the file to read.
    !! @return The model read from the file.
    function f_model_read(filepath) result(model)

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
    subroutine s_write_stl(filepath, model)

        character(LEN=*), intent(in) :: filepath
        type(t_model), intent(in) :: model

        integer :: i, j, iunit, iostat

        character(kind=c_char, len=80), parameter :: header = "Model file written by MFC."
        integer(kind=c_int32_t) :: nTriangles
        real(kind=c_float) :: normal(3), v(3)
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
    subroutine s_write_obj(filepath, model)

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
    subroutine s_model_write(filepath, model)

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

    function f_read_line(iunit, line) result(bIsLine)

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

    subroutine s_skip_ignored_lines(iunit)

        integer, intent(in) :: iunit

        character(80) :: line

        if (f_read_line(iunit, line)) then
            backspace (iunit)
        end if

    end subroutine s_skip_ignored_lines

    !> This procedure, recursively, finds whether a point is inside an octree.
    !! @param model    Model to search in.
    !! @param point    Point to test.
    !! @param spacing  Space around the point to search in (grid spacing).
    !! @param spc      Number of samples per cell.
    !! @return True if the point is inside the octree, false otherwise.
    function f_model_is_inside(model, point, spacing, spc) result(fraction)

        type(t_model), intent(in) :: model
        t_vec3, intent(in) :: point
        t_vec3, intent(in) :: spacing
        integer, intent(in) :: spc

        real(kind(0d0)) :: fraction

        type(t_ray) :: ray
        integer :: i, j, nInOrOut, nHits

        real(kind(0d0)), dimension(1:spc, 1:3) :: ray_origins, ray_dirs

        do i = 1, spc
            call random_number(ray_origins(i, :))
            ray_origins(i, :) = point + (ray_origins(i, :) - 0.5)*spacing(:)

            call random_number(ray_dirs(i, :))
            ray_dirs(i, :) = ray_dirs(i, :) - 0.5
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

    ! From https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
    !> This procedure checks if a ray intersects a triangle.
    !! @param ray      Ray.
    !! @param triangle Triangle.
    !! @return         True if the ray intersects the triangle, false otherwise.
    function f_intersects_triangle(ray, triangle) result(intersects)

        type(t_ray), intent(in) :: ray
        type(t_triangle), intent(in) :: triangle

        logical :: intersects

        real(kind(0d0)) :: v0v1(3), v0v2(3), N(3), P(3), C(3), edge(3), vp(3)
        real(kind(0d0)) :: area2, d, t, NdotRayDirection

        intersects = .false.

        N = triangle%n
        area2 = sqrt(sum(N(:)*N(:)))

        NdotRayDirection = sum(N(:)*ray%d(:))

        if (abs(NdotRayDirection) < 0.0000001) then
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

end module m_model
