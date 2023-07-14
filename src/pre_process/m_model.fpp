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
    subroutine s_read_stl_binary(filepath, model)

        character(LEN=*), intent(IN)  :: filepath
        type(t_model),    intent(OUT) :: model

        integer :: i, j, iunit, iostat

        character(kind=c_char, len=80) :: header
        integer  (kind=c_int32_t)      :: nTriangles

        real   (kind=c_float)          :: normal(3), v(3, 3)
        integer(kind=c_int16_t)        :: attribute

        open(newunit=iunit,      file=filepath, action='READ', &
             form='UNFORMATTED', status='OLD',  iostat=iostat, &
             access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open Binary STL file ", filepath
            
            call s_mpi_abort()
        end if

        read(iunit, iostat=iostat) header, nTriangles

        if (iostat /= 0) then
            print *, "Error: could not read header from Binary STL file ", filepath

            call s_mpi_abort()
        end if

        model%ntrs = nTriangles
        
        allocate(model%trs(model%ntrs))

        do i = 1, model%ntrs
            read(iunit) normal(:), v(1,:), v(2,:), v(3,:), attribute

            model%trs(i)%v = v
            model%trs(i)%n = normal
        end do

        close(iunit)

    end subroutine s_read_stl_binary

    !> This procedure reads an ASCII STL file.
    subroutine s_read_stl_ascii(filepath, model)

        character(LEN=*), intent(IN)  :: filepath
        type(t_model),    intent(OUT) :: model

        integer :: i, j, iunit, iostat

        character(80) :: line

        open(newunit=iunit,    file=filepath, action='READ', &
             form='FORMATTED', status='OLD',  iostat=iostat, &
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

        allocate(model%trs(model%ntrs))

        rewind(iunit)

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
            read(line(13:), *) model%trs(i)%n

            call s_skip_ignored_lines(iunit)
            read(iunit, '(A)') line

            do j = 1, 3
                if (.not. f_read_line(iunit, line)) exit

                if (line(1:6) /= "vertex") then
                    print *, "Error: expected vertex in STL file ", filepath

                    call s_mpi_abort()
                end if

                call s_skip_ignored_lines(iunit)
                read(line(7:), *) model%trs(i)%v(j,:)
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
    subroutine s_read_stl(filepath, model)

        character(LEN=*), intent(IN)  :: filepath
        type(t_model),    intent(OUT) :: model

        integer :: iunit, iostat

        character(80) :: line

        open(newunit=iunit,    file=filepath, action='READ', &
             form='FORMATTED', status='OLD',  iostat=iostat, &
             access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath
            
            call s_mpi_abort()
        end if

        read(iunit, '(A)') line
        
        close(iunit)

        if (line(1:5) == "solid") then
            call s_read_stl_ascii(filepath, model)
        else
            call s_read_stl_binary(filepath, model)
        end if
    
    end subroutine

    !> This procedure reads an OBJ file.
    subroutine s_read_obj(filepath, model)

        character(LEN=*), intent(IN)  :: filepath
        type(t_model),    intent(OUT) :: model

        integer :: i, j, k, l, iunit, iostat, nVertices

        t_vec3, allocatable :: vertices(:,:)
        
        character(80) :: line

        open(newunit=iunit,    file=filepath, action='READ', &
             form='FORMATTED', status='OLD',  iostat=iostat, &
             access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open model file ", filepath
            
            call s_mpi_abort()
        end if

        nVertices  = 0
        model%ntrs = 0
        do
            if (.not. f_read_line(iunit, line)) exit        

            select case(line(1:2))
            case ("v ")
                nVertices = nVertices + 1
            case ("f ")
                model%ntrs = model%ntrs + 1
            end select
        end do

        rewind(iunit)

        allocate(vertices(nVertices, 1:3))
        allocate(model%trs(model%ntrs))

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
                read(line(3:), *) vertices(i,:)
                i = i + 1
            case ("f ")
                read(line(3:), *) k, l, j
                model%trs(j)%v(1,:) = vertices(k,:)
                model%trs(j)%v(2,:) = vertices(l,:)
                model%trs(j)%v(3,:) = vertices(j,:)
                j = j + 1
            case default
                print *, "Error: unknown line type in OBJ file ", filepath
                print *, "Line: ", line

                call s_mpi_abort()
            end select
        end do

        deallocate(vertices)

        close(iunit)
 
    end subroutine

    !> This procedure reads a mesh from a file.
    !! @param filepath Path to the file to read.
    !! @return The model read from the file.
    function f_model_read(filepath) result(model)
    
        character(LEN=*), intent(IN) :: filepath

        type(t_model) :: model

        select case (filepath(len(trim(filepath))-3:len(trim(filepath))))
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
    subroutine s_write_stl(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model),    intent(IN) :: model

        integer :: i, j, iunit, iostat

        character(kind=c_char, len=80), parameter :: header = "Model file written by MFC."
        integer  (kind=c_int32_t) :: nTriangles
        real     (kind=c_float)   :: normal(3), v(3)
        integer  (kind=c_int16_t) :: attribute

        open(newunit=iunit,      file=filepath, action='WRITE', &
             form='UNFORMATTED', iostat=iostat, access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath
            
            call s_mpi_abort()
        end if

        nTriangles = model%ntrs
        write(iunit, iostat=iostat) header, nTriangles

        if (iostat /= 0) then
            print *, "Error: could not write header to STL file ", filepath

            call s_mpi_abort()
        end if

        do i = 1, model%ntrs
            normal = model%trs(i)%n
            write(iunit) normal

            do j = 1, 3
                v = model%trs(i)%v(j,:)
                write(iunit) v(:)
            end do

            attribute = 0
            write(iunit) attribute
        end do

        close(iunit)

    end subroutine s_write_stl

    !> This procedure writes an OBJ file.
    subroutine s_write_obj(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model),    intent(IN) :: model

        integer :: iunit, iostat

        integer :: i, j

        open(newunit=iunit,    file=filepath, action='WRITE', &
             form='FORMATTED', iostat=iostat, access='STREAM')
    
        if (iostat /= 0) then
            print *, "Error: could not open OBJ file ", filepath
            
            call s_mpi_abort()
        end if

        write(iunit, '(A)') "# Model file written by MFC."

        do i = 1, model%ntrs
            do j = 1, 3
                write(iunit, '(A, " ", (f30.20), " ", (f30.20), " ", (f30.20))') &
                    "v", model%trs(i)%v(j,1), model%trs(i)%v(j,2), model%trs(i)%v(j,3)
            end do

            write(iunit, '(A, " ", I0, " ", I0, " ", I0)') &
                "f", i*3-2, i*3-1, i*3
        end do

        close(iunit)

    end subroutine s_write_obj

    !> This procedure writes a binary STL file.
    !! @param filepath  Path to the file to write.
    !! @param triangles Triangles to write.
    subroutine s_model_write(filepath, model)
    
        character(LEN=*), intent(IN) :: filepath
        type(t_model),    intent(IN) :: model

        select case (filepath(len(trim(filepath))-3:len(trim(filepath))))
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

        type(t_model), intent(INOUT) :: model

        deallocate(model%trs)

    end subroutine s_model_free

    function f_read_line(iunit, line) result(bIsLine)

        integer,       intent(IN)  :: iunit
        character(80), intent(OUT) :: line
        logical :: bIsLine

        integer :: iostat

        bIsLine = .true.

        do
            read(iunit, '(A)', iostat=iostat) line

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

        integer, intent(IN) :: iunit

        character(80) :: line

        if (f_read_line(iunit, line)) then
            backspace(iunit)
        end if
        
    end subroutine s_skip_ignored_lines

end module m_model
