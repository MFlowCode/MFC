import rich

import os

import common, mfc

class MFCClean:
    def __init__(self, mfc: mfc.MFCState) -> None:
        self.mfc = mfc

    def clean_target(self, name: str, prepend="> ") -> None:
        if self.mfc.conf.get_target(name).fetch.method == "collection" or \
           self.mfc.args["recursive"]:
            sub_names = self.mfc.conf.get_dependency_names(name, recursive=self.mfc.args["recursive"])

            if len(sub_names) != 0:
                rich.print(f"{prepend}Package [bold blue]{name}[/bold blue] depends on {', '.join([ f'[bold blue]{x}[/bold blue]' for x in sub_names ])}:")

                for sub_name in sub_names:
                    self.clean_target(sub_name, prepend=prepend+"> ")

        rich.print(f"{prepend}Cleaning [bold blue]{name}[/bold blue]...")

        for clean_cmd in self.mfc.conf.get_target(name).clean:
            clean_cmd = self.mfc.build.string_replace(name, clean_cmd)

            src_path = self.mfc.build.get_source_path(name)

            if not self.mfc.lock.does_target_exist(name, self.mfc.args["mode"]):
                rich.print(f"{prepend}> This target was never built.")
                continue

            if os.path.exists(src_path):
                common.execute_shell_command(f'cd "{self.mfc.build.get_source_path(name)}" && {clean_cmd} > /dev/null 2>&1')
                # Set "bCleaned" flag
                self.mfc.lock.get_target(name, self.mfc.build. self.mfc.conf.get_desired_target_mode_name(name)).metadata.bCleaned=True
                self.mfc.lock.save()

    def run(self) -> None:
        for target_name in self.mfc.args["targets"]:
            self.clean_target(target_name)
