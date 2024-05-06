#!/usr/bin/env python3

import os, argparse

def main():
    parser = argparse.ArgumentParser(
        prog='indenter.py',
        description='Adjust indentation of OpenACC directives in a Fortran file')
    parser.add_argument('filepath', metavar='input_file', type=str, help='File to format')
    args = vars(parser.parse_args())

    filepath = args['filepath']

    temp_filepath = f"{filepath}.new"
    adjust_indentation(filepath, temp_filepath)
    os.replace(temp_filepath, filepath)

BLOCK_STARTERS  = ('if', 'do', '#:if', '#:else', "#ifdef", "#else")
BLOCK_ENDERS    = ('end', 'contains', 'else', '#:end', '#:else', '#else', '#endif')
LOOP_DIRECTIVES = ('!$acc loop', '!$acc parallel loop')
INDENTERS       = ('!DIR', '!$acc')

# pylint: disable=too-many-branches
def adjust_indentation(input_file, output_file):
    max_empty_lines=4
    indent_len=4

    with open(input_file, 'r') as file_in, open(output_file, 'w') as file_out:
        lines = file_in.readlines()

        # this makes sure !$acc lines that have line continuations get indented at proper level
        # pylint: disable=too-many-nested-blocks
        for _ in range(10):
            # loop through file
            # pylint: disable=consider-using-enumerate
            for i in range(len(lines)):
                if lines[i].lstrip().startswith(INDENTERS) and i + 1 < len(lines):
                    j = i + 1
                    empty_lines = 0
                    # look down to see how to indent a line
                    while j < len(lines) and empty_lines < max_empty_lines:
                        # if the following line starts with [end, else, contains], skip to looking up
                        if lines[j].lstrip().startswith(BLOCK_ENDERS):
                            empty_lines = max_empty_lines
                        # skip empty lines
                        elif lines[j].strip() == '':
                            empty_lines += 1
                        # indent acc lines
                        elif not lines[j].lstrip().startswith(INDENTERS):
                            indent = len(lines[j]) - len(lines[j].lstrip())
                            lines[i] = ' ' * indent + lines[i].lstrip()
                            break
                        j += 1
                    # if looking down just finds empty lines, start looking up for indentation level
                    if empty_lines == max_empty_lines:
                        k = i - 1
                        while k >= 0:
                            # if line above is not empty
                            if lines[k].strip() != '':
                                # if line 2 above ends with line continuation, indent at that level
                                if lines[k-1].strip().endswith('&'):
                                    indent = len(lines[k-1]) - len(lines[k-1].lstrip())
                                # if line above starts a loop or branch, indent
                                elif lines[k].lstrip().startswith(BLOCK_STARTERS):
                                    indent = indent_len + (len(lines[k]) - len(lines[k].lstrip()))
                                # else indent at level of line above
                                else:
                                    indent = len(lines[k]) - len(lines[k].lstrip())
                                lines[i] = ' ' * indent + lines[i].lstrip()
                                break
                            k -= 1

        # remove empty lines following an acc loop directive
        i = 0
        while i < len(lines):
            if lines[i].lstrip().startswith(LOOP_DIRECTIVES) and \
               i+1 < len(lines) and lines[i+1].strip() == '':
                file_out.write(lines[i])
                i += 2
            else:
                file_out.write(lines[i])
                i += 1

if __name__ == "__main__":
    main()
