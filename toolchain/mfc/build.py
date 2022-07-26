from mfc.util.printer import cons

import mfc.util.common as common

import os
import typing
import dataclasses


@dataclasses.dataclass
class MFCTarget:
    name:         str
    flags:        typing.List[str]
    requires:     typing.List[str]
    isDependency: bool


TARGETS: typing.List[MFCTarget] = [
    MFCTarget(name='fftw3', flags=['-DMFC_BUILD_FFTW3=ON'],
              isDependency=True, requires=[]),
    MFCTarget(name='hdf5', flags=['-DMFC_BUILD_HDF5=ON'],
              isDependency=True, requires=[]),
    MFCTarget(name='silo', flags=['-DMFC_BUILD_SILO=ON'],
              isDependency=True, requires=["hdf5"]),
    MFCTarget(name='pre_process', flags=['-DMFC_BUILD_PRE_PROCESS=ON'],
              isDependency=False, requires=[]),
    MFCTarget(name='simulation', flags=['-DMFC_BUILD_SIMULATION=ON'],
              isDependency=False, requires=["fftw3"]),
    MFCTarget(name='post_process', flags=['-DMFC_BUILD_POST_PROCESS=ON'],
              isDependency=False, requires=['fftw3', 'silo'])
]


def get_mfc_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS if not target.isDependency ]


def get_dependencies_names() -> typing.List[str]:
    return [ target.name for target in TARGETS if target.isDependency ]


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


def is_target_configured(target: MFCTarget) -> bool:
    build_dirpath = get_build_dirpath(target)
    return os.path.isdir(build_dirpath)


def clean_target(mfc, name: str):
    cons.print(f"Cleaning [bold magenta]{name}[/bold magenta]:")
    cons.indent()

    target = get_target(name)
    mode   = mfc.user.get_mode(mfc.args["mode"])

    build_dirpath = get_build_dirpath(target)

    if not os.path.isdir(build_dirpath):
        cons.print("Target not configured. Nothing to clean.")
        cons.unindent()
        return

    clean = f"cmake --build . --target clean --config {mode.type}"

    cons.print()
    cons.print(f"$ {clean}", highlight=False)
    cons.print()
    common.system(f"cd \"{build_dirpath}\" && {clean}", exception_text=f"Failed to clean the [bold magenta]{name}[/bold magenta] target.")

    cons.unindent()


def build_target(mfc, name: str, history: typing.List[str] = None):
    cons.print(f"[bold]Building [magenta]{name}[/magenta]:[/bold]")
    cons.indent()

    if history is None:
        history = []

    if name in history:
        cons.print("Already built, skipping...")
        cons.unindent()
        return

    history.append(name)

    target = get_target(name)
    mode   = mfc.user.get_mode(mfc.args["mode"])

    if target.isDependency and mfc.args[f"no_{target.name}"]:
        cons.print(f"--no-{target.name} given, skipping...")
        cons.unindent()
        return

    build_dirpath   = get_build_dirpath(target)
    cmake_dirpath   = get_cmake_dirpath(target)
    install_dirpath = get_install_dirpath()

    flags: list = target.flags.copy() + mode.flags + [
        f"-Wno-dev",
        f"-DCMAKE_BUILD_TYPE={mode.type}",
        f"-DCMAKE_PREFIX_PATH=\"{install_dirpath}\"",
        f"-DCMAKE_INSTALL_PREFIX=\"{install_dirpath}\"",
    ]

    configure = f"cd \"{build_dirpath}\" && cmake {' '.join(flags)} \"{cmake_dirpath}\""
    build     = f"cd \"{build_dirpath}\" && cmake --build . -j {mfc.args['jobs']} --target {name} --config {mode.type}"
    install   = f"cd \"{build_dirpath}\" && cmake --install ."

    # Only configure the first time
    if not is_target_configured(target):
        for dependency_name in target.requires:
            build_target(mfc, dependency_name, history)

        cons.print(no_indent=True)
        cons.print(f"$ {configure}")
        cons.print(no_indent=True)
        common.create_directory(build_dirpath)

        if common.system(configure, no_exception=True) != 0:
            common.delete_directory(build_dirpath)

            raise common.MFCException(f"Failed to configure the [bold magenta]{name}[/bold magenta] target.")

    cons.print(no_indent=True)
    cons.print(f"$ {build}")
    cons.print(no_indent=True)
    common.system(build, exception_text=f"Failed to build the [bold magenta]{name}[/bold magenta] target.")
    cons.print(no_indent=True)
    cons.print(f"$ {install}")
    cons.print(no_indent=True)
    common.system(install, exception_text=f"Failed to install the [bold magenta]{name}[/bold magenta] target.")
    cons.print(no_indent=True)

    cons.unindent()


def build(mfc):
    cons.print("[bold]Build:[/bold]")
    cons.indent()

    for target_name in mfc.args["targets"]:
        build_target(mfc, target_name)

    cons.unindent()
