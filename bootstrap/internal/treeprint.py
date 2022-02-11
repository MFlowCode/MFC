from dataclasses import dataclass
import internal.common as common

import colorama
import dataclasses

@dataclasses.dataclass
class Indent:
    color: str
    s:     str

class TreePrinter:
    def __init__(self):
        self.indents = [Indent("", "├─ ")]
    
    def get_depth(self) -> int:
        return len(self.indents)

    def indent(self, color=''):
        if len(self.indents) == 0:
            self.indents.append(Indent(color, '├─ '))
        else:
            self.indents[-1].s = '|  '
            self.indents.append(Indent(color, '├─ '))
        
    def unindent(self):
        self.indents.pop()

        if len(self.indents) != 0:
            self.indents[-1].s = '├─ '
    
    def get_format(self, s, end='\n'):
        indent_str = ""
        for indent in self.indents:
            indent_str += f"{indent.color}{indent.s}{colorama.Style.RESET_ALL}"
        
        return indent_str + s + end

    def print(self, s, end='\n'):
        print(self.get_format(s, end=end), end='')

    def print_progress(self, s, i, e):
        length  = 20
        nFilled = int((i*length)/e)
        nBlank  = length - nFilled

        end = '\r' if i != e else '\n'

        common.clear_line()
        self.print(f'{s} [{nFilled*"#"}{nBlank*" "}] {i}/{e} - {int(i/e*100)}%', end=end)

