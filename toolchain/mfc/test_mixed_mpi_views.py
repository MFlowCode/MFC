"""Exercise production MPI file views and I/O calls without requiring REAL(2)."""

import os
import re
import shutil
import struct
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

# The production initializers only associate these buffers, so opaque two-byte
# cells reproduce half storage on compilers that cannot do half arithmetic.
GLOBALS = """
module test_globals
    use mpi
    use iso_c_binding, only: c_int8_t, c_int16_t
    implicit none
    integer, parameter :: mpi_p=MPI_REAL, mpi_io_p=MPI_BYTE, mpi_io_type=2
    type, bind(c) :: cell_t
        integer(c_int8_t) :: bits(2)
    end type
    type scalar_field
        type(cell_t), pointer :: sf(:,:,:) => null()
    end type
    type pres_field
        type(cell_t), pointer :: sf(:,:,:,:,:) => null()
    end type
    type integer_field
        integer, pointer :: sf(:,:,:) => null()
    end type
    type mpi_io_ib_var
        integer :: view
        type(integer_field) :: var
    end type
    type mpi_io_var
        integer, allocatable :: view(:)
        type(scalar_field), allocatable :: var(:)
    end type
    type(mpi_io_var) :: MPI_IO_DATA
    integer :: m,n,p,m_glb,n_glb,p_glb,num_dims,sys_size,nb,nnode,start_idx(3)
    logical :: qbmm,polytropic
contains
    subroutine s_mpi_abort(message)
        character(len=*), intent(in) :: message
        integer :: ierr
        print *, message
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end subroutine
end module
"""

DRIVER = """
program mixed_mpi_views
    use test_globals
    use production_views
    implicit none
    integer :: rank, ranks, ierr, ifile, data_size, mpi_info_int, i, x, y, writer, reader
    integer :: status(MPI_STATUS_SIZE), errors, total_errors, guard_errors, total_guard_errors, code
    integer(kind=MPI_OFFSET_KIND) :: disp
    type(scalar_field) :: q(3)
    type(cell_t), target :: backing(18,3)
    integer(c_int8_t) :: expected(2)
    character(len=32) :: filename

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ranks, ierr)
    if (ranks /= 2) call s_mpi_abort('This test requires two ranks.')
    num_dims=2; sys_size=3; nb=1; nnode=2
    m=2; n=2; p=0; m_glb=2; n_glb=5; p_glb=0
    start_idx=[0,rank*3,0]
    qbmm=.false.; polytropic=.true.
    data_size=9
    mpi_info_int=MPI_INFO_NULL
    allocate(MPI_IO_DATA%var(sys_size),MPI_IO_DATA%view(sys_size))
    do i=1,sys_size
        ! Nine guard cells absorb and reveal an erroneous four-byte read.
        q(i)%sf(0:m,0:n,0:p) => backing(1:data_size,i)
    end do
    call s_initialize_mpi_data(q)
    errors=0; guard_errors=0
    do writer=0,__WRITERS_COUNT__
        filename='serial.dat'
        if (writer > 0) then
            write(filename,'(A,I0,A)') 'writer',writer,'.dat'
            do i=1,sys_size
                do y=0,n
                    do x=0,m
                        code=i*1000+(start_idx(2)+y)*3+x+1
                        q(i)%sf(x,y,0)%bits=transfer(int(code,c_int16_t),expected)
                    end do
                end do
            end do
            call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_CREATE+MPI_MODE_RDWR,mpi_info_int,ifile,ierr)
            call check_mpi()
            call MPI_FILE_SET_SIZE(ifile,0_MPI_OFFSET_KIND,ierr)
            call check_mpi()
            do i=1,sys_size
                disp=int(i-1,MPI_OFFSET_KIND)*18*2
                select case(writer)
__WRITERS__
                end select
            end do
            call MPI_FILE_CLOSE(ifile,ierr)
            call check_mpi()
        end if
        ! Read both the independent serial reference and each parallel writer.
        do reader=1,__READERS_COUNT__
            do i=1,sys_size
                do x=1,18
                    backing(x,i)%bits=-77_c_int8_t
                end do
            end do
            call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,mpi_info_int,ifile,ierr)
            call check_mpi()
            do i=1,sys_size
                disp=int(i-1,MPI_OFFSET_KIND)*18*2
                select case(reader)
__READERS__
                end select
                do y=0,n
                    do x=0,m
                        code=i*1000+(start_idx(2)+y)*3+x+1
                        expected=transfer(int(code,c_int16_t),expected)
                        errors=errors+count(q(i)%sf(x,y,0)%bits /= expected)
                    end do
                end do
                do x=data_size+1,18
                    errors=errors+count(backing(x,i)%bits /= -77_c_int8_t)
                    guard_errors=guard_errors+count(backing(x,i)%bits /= -77_c_int8_t)
                end do
            end do
            call MPI_FILE_CLOSE(ifile,ierr)
            call check_mpi()
        end do
    end do
    call MPI_ALLREDUCE(errors,total_errors,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(guard_errors,total_guard_errors,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (rank == 0) print *, 'mismatched data or guard bytes:',total_errors
    if (rank == 0) print *, 'modified guard bytes:',total_guard_errors
    do i=1,sys_size
        call MPI_TYPE_FREE(MPI_IO_DATA%view(i),ierr)
    end do
    call MPI_FINALIZE(ierr)
    if (total_errors /= 0) stop 1
contains
    subroutine check_mpi(operation)
        character(len=*), optional, intent(in) :: operation
        character(len=MPI_MAX_ERROR_STRING) :: message
        integer :: abort_ierr, length
        if (ierr == MPI_SUCCESS) return
        call MPI_ERROR_STRING(ierr,message,length,abort_ierr)
        print *, 'MPI error on rank',rank,': ',message(:length)
        if (present(operation)) print *, 'operation: ',operation
        call MPI_ABORT(MPI_COMM_WORLD,ierr,abort_ierr)
    end subroutine
end program
"""


def production_initializers():
    """Compile the actual initializer bodies, including their helper if present."""
    source = (ROOT / "src/common/m_mpi_common.fpp").read_text()
    bodies = []
    for name in ("s_initialize_mpi_data", "s_create_mpi_field_subarray"):
        match = re.search(r"^    (?:impure )?subroutine " + name + r"\(.*?^    end subroutine " + name + r"\s*$", source, re.M | re.S)
        if match is None:
            assert name != "s_initialize_mpi_data", "production MPI initializer not found"
            continue
        body = match.group()
        # These routines' only preprocessor condition is MFC_MPI. Enable it
        # without imposing compiler-specific preprocessing flags on the fixture.
        body = re.sub(r"^#(?:ifdef MFC_MPI|endif)\s*$", "", body, flags=re.M)
        assert not re.search(r"^#", body, re.M), "unexpected initializer preprocessing"
        bodies.append(body)
    return "module production_views\nuse test_globals\nimplicit none\ncontains\n" + "\n".join(bodies) + "\nend module\n"


def production_io_cases(paths, operation):
    """Use literal view/transfer calls so wrong production count/type is tested."""
    cases = []
    for path in paths:
        source = (ROOT / path).read_text()
        pairs = re.findall(
            r"(call MPI_FILE_SET_VIEW\(ifile, disp, [^\n]+MPI_IO_DATA%view\(i\)[^\n]+\))\s*\n"
            r"\s*(call MPI_FILE_" + operation + r"(?:_ALL)?\(ifile, MPI_IO_DATA%var\(i\)%sf(?:\(:, ?:, ?:\))?,[^\n]+\))",
            source,
        )
        assert pairs, f"no production {operation} call found in {path}"
        for view, transfer in dict.fromkeys(pairs):
            cases.append(f"case({len(cases) + 1})\n{view}\ncall check_mpi('SET_VIEW')\n{transfer}\ncall check_mpi('{operation}')")
    return "\n".join(cases), len(cases)


def test_mixed_parallel_file_views(tmp_path):
    if os.environ.get("MFC_SKIP_COMPILER_TESTS", "") == "1":
        pytest.skip("MFC_SKIP_COMPILER_TESTS=1 — skipping tests that compile Fortran")
    compiler = shutil.which("mpifort") or shutil.which("mpif90")
    launcher = shutil.which("mpiexec") or shutil.which("mpirun")
    if compiler is None or launcher is None:
        pytest.skip("an MPI Fortran compiler and MPI launcher are required")

    writers = ["src/pre_process/m_data_output.fpp", "src/simulation/m_data_output.fpp"]
    readers = ["src/pre_process/m_start_up.fpp", "src/simulation/m_start_up.fpp", "src/post_process/m_data_input.f90"]
    write_calls, writer_count = production_io_cases(writers, "WRITE")
    read_calls, reader_count = production_io_cases(readers, "READ")
    driver = DRIVER.replace("__WRITERS__", write_calls).replace("__READERS__", read_calls)
    driver = driver.replace("__READERS_COUNT__", str(reader_count)).replace("__WRITERS_COUNT__", str(writer_count))
    sources = [tmp_path / name for name in ("globals.f90", "views.f90", "driver.f90")]
    for path, source in zip(sources, (GLOBALS, production_initializers(), driver)):
        path.write_text(source)
    executable = tmp_path / "mpi_views"
    compiled = subprocess.run([compiler, *map(str, sources), "-o", str(executable)], cwd=tmp_path, capture_output=True, text=True, timeout=60, check=False)
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr

    # Three fields on a 3x6 grid, split into two 3x3 (odd-sized) blocks.
    expected = b"".join(struct.pack("=h", field * 1000 + cell) for field in range(1, 4) for cell in range(1, 19))
    (tmp_path / "serial.dat").write_bytes(expected)
    env = os.environ.copy()
    # Open MPI otherwise rejects container CI users and two ranks on one slot.
    env.update(OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1", OMPI_MCA_rmaps_base_oversubscribe="1")
    result = subprocess.run([launcher, "-n", "2", str(executable)], cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    for writer in range(1, writer_count + 1):
        assert (tmp_path / f"writer{writer}.dat").read_bytes() == expected, f"writer {writer} disagrees with the serial reference"
