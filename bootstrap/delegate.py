#!/usr/bin/env python3

import colorama
import traceback

import internal.common    as common
import internal.bootstrap as bootstrap

def main():
    common.execute_shell_command_safe("git submodule update --init --recursive", no_exception=True)

    try:
        common.colorama.init()

        bootstrap.Bootstrap()
    except common.MFCException as exc:
        print(traceback.format_exc())
        print(f"{colorama.Fore.RED}|--> {str(exc)}{colorama.Style.RESET_ALL}")
        exit(1)


if __name__ == "__main__":
    main()

