import os, dataclasses

import conf, common

@dataclasses.dataclass
class LockTargetMetadata:
    mode:     str
    bCleaned: bool

    def __init__(self, data: dict):
        self.mode     = data["mode"]
        self.bCleaned = data["bCleaned"]

@dataclasses.dataclass
class LockTargetHolder:
    target:   conf.Target
    metadata: LockTargetMetadata

    def __init__(self, data: dict):
        self.target   = conf.Target(data["target"])
        self.metadata = LockTargetMetadata(data["metadata"])

@dataclasses.dataclass
class MFCLock:
    targets: list
    mode:    str

    def __init__(self, mfc):
        default_mode: str = mfc.user.modes[0].name

        if not os.path.exists(common.MFC_LOCK_FILEPATH):
            common.create_file(common.MFC_LOCK_FILEPATH)
            common.file_dump_yaml(common.MFC_LOCK_FILEPATH, {
                "mode": default_mode,
                "targets": []
            })

        self.data: dict = common.file_load_yaml(common.MFC_LOCK_FILEPATH)

        self.targets = []
        self.mode    = self.data.get("mode", default_mode)

        for t in self.data["targets"]:
            self.add_target(LockTargetHolder(t))

    def save(self):
        self.flush()
        common.file_dump_yaml(common.MFC_LOCK_FILEPATH, self.data)

    def flush(self):
        self.data = dataclasses.asdict(self)

    def add_target(self, target: LockTargetHolder):
        self.targets.append(target)
        self.flush()

    def get_target_matches(self, name: str, restrict_cc: str = None) -> list:
        def peek_filter(e: dict):
            if e.target.name != name:
                return False

            if restrict_cc is None:
                return True

            return e.metadata.mode == restrict_cc

        return list(filter(peek_filter, self.targets))

    def does_target_exist(self, name: str, restrict_cc: str = None) -> bool:
        return len(self.get_target_matches(name, restrict_cc)) > 0

    def does_unique_target_exist(self, name: str, restrict_cc: str = None) -> bool:
        return len(self.get_target_matches(name, restrict_cc)) == 1

    def get_target(self, name: str, restrict_cc: str = None) -> LockTargetHolder:
        matches = self.get_target_matches(name, restrict_cc)

        if len(matches) == 0:
            raise common.MFCException(f'Failed to retrieve dependency "{name}" with restrict_cc="{restrict_cc}".')

        if len(matches) > 1:
            raise common.MFCException(f'More than one dependency to choose from for "{name}" with restrict_cc="{restrict_cc}".')

        return matches[0]
