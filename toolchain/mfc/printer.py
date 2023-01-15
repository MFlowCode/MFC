import typing

import rich, rich.console


class MFCPrinter:
    def __init__(self):
        self.stack = []
        self.raw   = rich.console.Console()

    def reset(self):
        self.stack = []

    def indent(self, msg: str = None):
        msg = msg if msg is not None else "  "
        
        self.stack.append(msg)

    def unindent(self, times: int = None):
        if times == None:
            times = 1
        
        for _ in range(times):
            self.stack.pop()

    def print(self, msg: typing.Any = None, no_indent: bool = False, *args, **kwargs):
        if msg is None:
            msg = ""
        
        if no_indent:
            self.raw.print(str(msg), *args, **kwargs)
        else:
            print_s, lines = "", str(msg).split('\n', maxsplit=-1)
            for i, s in enumerate(lines):
                newline = '\n' if (i != len(lines)-1) else ''
                print_s += f"{''.join(self.stack)}{s}{newline}"

            self.raw.print(print_s, *args, **kwargs)

    def print_exception(self):
        self.raw.print_exception()


cons = MFCPrinter()
