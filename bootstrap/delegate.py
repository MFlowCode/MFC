#!/usr/bin/env python3

import sys
import rich
import signal
import traceback

import internal.common    as common
import internal.bootstrap as bootstrap

def main():
    try:
        bootstrap.Bootstrap()
    except common.MFCException as exc:
        rich.print(f"[red]> {str(exc)}[/red]")
        exit(1)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt as exc:
        exit(0)
