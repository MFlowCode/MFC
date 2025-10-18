import dataclasses, typing, sys, os, re, math
from datetime import datetime
from pathlib import Path

from ..       import common
from ..run    import input
from ..build  import get_configured_targets
from ..state  import CFG




# This class maps to the data contained in one file in D/
@dataclasses.dataclass(repr=False)
class PackEntry:
    filepath: str
    doubles:  typing.List[float]

    def __repr__(self) -> str:
        return f"{self.filepath} {' '.join([ str(d) for d in self.doubles ])}"


# This class maps to the data contained in the entirety of D/: it is tush a list
# of PackEntry classes.
class Pack:
    entries: typing.Dict[str, PackEntry]

    def __init__(self, entries: typing.List[PackEntry] = None):
        self.entries = {}

        for entry in entries or []:
            self.set(entry)

    def find(self, filepath: str) -> PackEntry:
        return self.entries.get(filepath, None)

    def set(self, entry: PackEntry):
        self.entries[entry.filepath] = entry

    def remove(self, filepath_or_entry: typing.Union[str, PackEntry]):
        if isinstance(filepath_or_entry, str):
            del self.entries[filepath_or_entry]
        else:
            del self.entries[filepath_or_entry.filepath]

    def save(self, filepath: str):
        if filepath.endswith(".py"):
            filepath = os.path.dirname(filepath)

        if os.path.isdir(filepath):
            filepath = os.path.join(filepath, "pack.txt")

        if not filepath.endswith(".txt"):
            filepath += ".txt"

        common.file_write(filepath, '\n'.join([ str(e) for e in sorted(self.entries.values(), key=lambda x: x.filepath) ]))

        metadata = f"""\
This file was created on {str(datetime.now())}.

mfc.sh:

    Invocation: {' '.join(sys.argv[1:])}
    Lock:       {CFG()}
    Git:        {common.generate_git_tagline()}

"""

        case = input.load(None, {}, {})
        for target in get_configured_targets(case):
            cfg = target.get_configuration_txt(case)

            if cfg is None:
                continue

            metadata += f"""\
{target.name}:

    {'    '.join(cfg.splitlines(keepends=True))}
"""

        metadata += f"""\
CPU:

    {'    '.join(common.get_cpuinfo().splitlines(keepends=True))}
"""

        common.file_write(f"{filepath.rstrip('.txt')}-metadata.txt", metadata)

    def has_NaNs(self) -> bool:
        for entry in self.entries.values():
            for double in entry.doubles:
                if math.isnan(double):
                    return True

        return False


def load(filepath: str) -> Pack:
    if not os.path.isfile(filepath):
        filepath = os.path.join(filepath, "pack.txt")

    entries: typing.List[PackEntry] = []

    for line in common.file_read(filepath).splitlines():
        if common.isspace(line):
            continue

        arr = line.split(' ')

        entries.append(PackEntry(
            filepath=arr[0],
            doubles=[ float(d) for d in arr[1:] ]
        ))

    return Pack(entries)


def compile(casepath: str) -> typing.Tuple[Pack, str]:
    entries = []

    case_dir = os.path.dirname(casepath) if os.path.isfile(casepath) else casepath
    D_dir    = os.path.join(case_dir, "D")

    for filepath in list(Path(D_dir).rglob("*.dat")):
        short_filepath = str(filepath).replace(f'{case_dir}', '')[1:].replace("\\", "/")
        content        = common.file_read(filepath)

        # Takes a string of numbers and returns them as a list of floats.
        def _extract_doubles(s: str) -> list:
            return [ float(e) for e in re.sub(r"[\n\t\s]+", " ", s).strip().split(' ') ]


        try:
            if "lag_bubble" in short_filepath:
                lines = content.splitlines()
                content = "\n".join(lines[1:])  # Skip the first line
                doubles = _extract_doubles(content)
            else:
                # Every line is <x> <y> <z> <value> (<y> and <z> are optional). So the
                # number of dimensions is the number of doubles in the first line minus 1.
                ndims   = len(_extract_doubles(content.split('\n', 1)[0])) - 1
                # We discard all <x> <y> <z> values and only keep the <value> ones.
                # This is in an effort to save on storage.
                doubles = _extract_doubles(content)[ndims::ndims+1]
        except ValueError:
            return None, f"Failed to interpret the content of [magenta]{filepath}[/magenta] as a list of floating point numbers."

        entries.append(PackEntry(short_filepath,doubles))

    return Pack(entries), None
