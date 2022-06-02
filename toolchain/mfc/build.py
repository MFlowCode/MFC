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
    isCollection: bool
    requires:     typing.List[str]


TARGETS: typing.List[MFCTarget] = [
    MFCTarget(
        name='fftw3',
        flags=['-DMFC_BUILD_FFTW3=ON'],
        isDependency=True,
        isCollection=False,
        requires=[]
    ), MFCTarget(
        name='silo',
        flags=['-DMFC_BUILD_SILO=ON'],
        isDependency=True,
        isCollection=False,
        requires=[]
    ), MFCTarget(
        name='pre_process',
        flags=['-DMFC_BUILD_PRE_PROCESS=ON'],
        isDependency=False,
        isCollection=False,
        requires=[]
    ), MFCTarget(
        name='simulation',
        flags=['-DMFC_BUILD_SIMULATION=ON'],
        isDependency=False,
        isCollection=False,
        requires=["fftw3"]
    ), MFCTarget(
        name='post_process',
        flags=['-DMFC_BUILD_POST_PROCESS=ON'],
        isDependency=False,
        isCollection=False,
        requires=['fftw3', 'silo']
    ), MFCTarget(
        name='mfc',
        flags=[],
        isDependency=False,
        isCollection=True,
        requires=['pre_process', 'simulation', 'post_process']
    ),
    MFCTarget(
        name="clean",
        flags=[],
        isDependency=False,
        isCollection=True,
        requires=[])
]


def get_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS ]


def get_target(name: str) -> MFCTarget:
    for target in TARGETS:
        if target.name == name:
            return target
    
    raise common.MFCException(f"Target '{name}' does not exist.")


# Get path to directory that will store the build files
def get_build_dirpath(target: MFCTarget) -> str:
    return f'{os.getcwd()}/build/{["mfc", "dependencies"][int(target.isDependency)]}/{target.name}'


# Get the directory that contains the target's CMakeLists.txt
def get_cmake_dirpath(target: MFCTarget) -> str:
    return f'{os.getcwd()}/{["", "toolchain/dependencies"][int(target.isDependency)]}'


def get_install_dirpath() -> str:
    return f'{os.getcwd()}/build/install'


def clean_target(mfc, name: str):
    target = get_target(name)

    if target.isCollection:
        for dependency_name in target.requires:
            clean_target(mfc, dependency_name)
        
        return

    build_dirpath = get_build_dirpath(target)
    clean_cmd     = f"cd '{build_dirpath}' && cmake --build . --target clean"

    rich.print(f" + {clean_cmd}")
    common.execute_shell_command(clean_cmd)


def build_target(mfc, name: str, history: typing.List[str] = None):
    if history is None:
        history = []

    if name in history:
        return
    
    history.append(name)

    target = get_target(name)

    for dependency_name in target.requires:
        build_target(mfc, dependency_name, history)

    if target.isCollection:
        return

    flags: list = target.flags.copy()

    if mfc.args["mode"] is not None:
        flags = flags + mfc.user.get_mode(mfc.args["mode"]).flags

    build_dirpath   = get_build_dirpath(target)
    cmake_dirpath   = get_cmake_dirpath(target)
    install_dirpath = get_install_dirpath()

    configure = f"cd '{build_dirpath}' && cmake -DCMAKE_INSTALL_PREFIX='{install_dirpath}' -DCMAKE_PREFIX_PATH='{install_dirpath}' {' '.join(flags)} '{cmake_dirpath}'"
    build     = f"cd '{build_dirpath}' && cmake --build . -j {mfc.args['jobs']} --target {name}"
    install   = f"cd '{build_dirpath}' && cmake --install ."

    # Only configure the first time
    if not os.path.exists(build_dirpath):
        os.mkdir(build_dirpath)
        rich.print(f" + {configure}")
        common.execute_shell_command(configure)

    rich.print(f" + {build}")
    common.execute_shell_command(build)

    rich.print(f" + {install}")
    common.execute_shell_command(install)
