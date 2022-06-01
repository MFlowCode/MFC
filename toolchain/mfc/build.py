import rich

import common

import os
import typing
import dataclasses


@dataclasses.dataclass
class MFCTarget:
    name:         str
    flags:        typing.List[str]
    isDependency: bool
    requires:     list[str]


TARGETS: typing.List[MFCTarget] = [
    MFCTarget(
        name='fftw3',
        flags=['-DMFC_BUILD_FFTW3=ON'],
        isDependency=True,
        requires=[]
    ), MFCTarget(
        name='hdf5',
        flags=['-DMFC_BUILD_HDF5=ON'],
        isDependency=True,
        requires=[]
    ), MFCTarget(
        name='silo',
        flags=['-DMFC_BUILD_SILO=ON'],
        isDependency=True,
        requires=['hdf5']
    ), MFCTarget(
        name='pre_process',
        flags=['-DMFC_BUILD_PRE_PROCESS=ON'],
        isDependency=False,
        requires=[]
    ), MFCTarget(
        name='simulation',
        flags=['-DMFC_BUILD_SIMULATION=ON'],
        isDependency=False,
        requires=["fftw3"]
    ), MFCTarget(
        name='post_process',
        flags=['-DMFC_BUILD_POST_PROCESS=ON'],
        isDependency=False,
        requires=['fftw3', 'hdf5', 'silo']
    ), MFCTarget(
        name='mfc',
        flags=[],
        isDependency=False,
        requires=['pre_process', 'simulation', 'post_process']
    ),
    MFCTarget(
        name="clean",
        flags=[],
        isDependency=False,
        requires=[])
]


def get_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS ]


def get_target(name: str) -> MFCTarget:
    for target in TARGETS:
        if target.name == name:
            return target
    
    raise common.MFCException(f"Target '{name}' does not exist.")


def clean(mfc):
    build_target(mfc, "clean")


def build_target(mfc, name: str):
    target = get_target(name)

    for dependency_name in target.requires:
        build_target(mfc, dependency_name)

    flags: list = target.flags.copy()

    if mfc.args["mode"] is not None:
        flags = flags + mfc.user.get_mode(mfc.args["mode"]).flags

    build_dirpath   = f'{os.getcwd()}/build/{["mfc", "dependencies"][int(target.isDependency)]}'
    cmake_dirpath   = f'{os.getcwd()}/{["", "toolchain/dependencies"][int(target.isDependency)]}'
    install_dirpath = f'{os.getcwd()}/build/install'

    cmake_target = name if name != 'mfc' else 'all'

    configure = f"cd '{build_dirpath}' && cmake -DCMAKE_INSTALL_PREFIX='{install_dirpath}' {' '.join(flags)} '{cmake_dirpath}'"
    build     = f"cd '{build_dirpath}' && cmake --build . -j {mfc.args['jobs']} --target {cmake_target}"
    install   = f"cd '{build_dirpath}' && cmake --install ."

    rich.print(f" + {configure}")
    common.execute_shell_command(configure)
    rich.print(f" + {build}")
    common.execute_shell_command(build)

    if name != "clean":
        rich.print(f" + {install}")
        common.execute_shell_command(install)
