import os, dataclasses

import conf, common

@dataclasses.dataclass
class LockTargetMetadata:
    bCleaned: bool
    mode: str

    def __init__(self, data: dict):
        self.bCleaned               = data.get("bCleaned", False)
        self.mode = data["mode"]

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

    def __init__(self):
        if not os.path.exists(common.MFC_LOCK_FILEPATH):
            common.create_file(common.MFC_LOCK_FILEPATH)
            common.file_dump_yaml(common.MFC_LOCK_FILEPATH, {"targets": []})

        self.data: dict = common.file_load_yaml(common.MFC_LOCK_FILEPATH)

        self.targets = []

        for t in self.data["targets"]:
            self.add_target(LockTargetHolder(t))

        self.flush()
    
    def save(self):
        common.file_dump_yaml(common.MFC_LOCK_FILEPATH, self.data)

    def flush(self):
        self.data = dataclasses.asdict(self)

    def add_target(self, target: LockTargetHolder):
        self.targets.append(target)
        self.flush()

    def was_target_built(self, name: str, restrict_cc: str = None):
        matches = self.get_target_matches(name, restrict_cc)

    def get_target_matches(self, name: str, restrict_cc: str = None):
        def peek_filter(e: dict):
            if e.target.name != name:
                return False

            if restrict_cc is None:
                return True

            return e.metadata.mode == restrict_cc

        return list(filter(peek_filter, self.targets))

    def does_target_exist(self, name: str, restrict_cc: str = None):
        return len(self.get_target_matches(name, restrict_cc)) > 0

    def does_unique_target_exist(self, name: str, restrict_cc: str = None):
        return len(self.get_target_matches(name, restrict_cc)) == 1

    def get_target(self, name: str, restrict_cc: str = None):
        matches = self.get_target_matches(name, restrict_cc)

        if len(matches) == 0:
            raise common.MFCException(f'Failed to retrieve dependency "{name}" with restrict_cc="{restrict_cc}".')

        if len(matches) > 1:
            raise common.MFCException(f'More than one dependency to choose from for "{name}" with restrict_cc="{restrict_cc}".')

        return matches[0]
