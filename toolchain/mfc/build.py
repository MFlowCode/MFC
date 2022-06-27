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
    return os.sep.join([os.getcwd(), "build", target.name])


# Get the directory that contains the target's CMakeLists.txt
def get_cmake_dirpath(target: MFCTarget) -> str:
    return os.sep.join([os.getcwd(), ["", "toolchain/dependencies"][int(target.isDependency)]])


def get_install_dirpath() -> str:
    return os.sep.join([os.getcwd(), "build", "install"])


def clean_target(mfc, name: str):
    target = get_target(name)
    mode   = mfc.user.get_mode(mfc.args["mode"])

    if target.isCollection:
        for dependency_name in target.requires:
            clean_target(mfc, dependency_name)
        
        return

    build_dirpath = get_build_dirpath(target)

    if not os.path.isdir(build_dirpath):
        return

    clean = f"cmake --build . --target clean --config {mode.type}"

    rich.print(f" + {clean}")
    common.system(f"cd \"{build_dirpath}\" && {clean}")


def build_target(mfc, name: str, history: typing.List[str] = None):
    if history is None:
        history = []

    if name in history:
        return
    
    history.append(name)

    target = get_target(name)
    mode   = mfc.user.get_mode(mfc.args["mode"])

    for dependency_name in target.requires:
        build_target(mfc, dependency_name, history)

    if target.isCollection:
        return

    build_dirpath   = get_build_dirpath(target)
    cmake_dirpath   = get_cmake_dirpath(target)
    install_dirpath = get_install_dirpath()

    flags: list = target.flags.copy() + mode.flags + [
        f"-Wno-dev"
        f"-DCMAKE_BUILD_TYPE={mode.type}",
        f"-DCMAKE_PREFIX_PATH=\"{install_dirpath}\"",
        f"-DCMAKE_INSTALL_PREFIX=\"{install_dirpath}\"",
    ]

    cd        = f"cd \"{build_dirpath}\""
    configure = f"cmake {' '.join(flags)} \"{cmake_dirpath}\""
    build     = f"cmake --build . -j {mfc.args['jobs']} --target {name} --config {mode.type}"
    install   = f"cmake --install ."

    # Only configure the first time
    if not os.path.exists(build_dirpath):
        rich.print(f" + {configure}")
        common.create_directory(build_dirpath)
        
        if common.system(f"{cd} && {configure}", no_exception=True) != 0:
            common.delete_directory(build_dirpath)

            raise common.MFCException(f"Failed to configure the {name} target.")
            
    rich.print(f" + {build}")
    common.system(f"{cd} && {build}")

    rich.print(f" + {install}")
    common.system(f"{cd} && {install}")
