#!/usr/bin/env python3

import rich
#import traceback

import user, conf, lock, args, test, clean, common, run

class MFCState:
    def __init__(self) -> None:
        from build import MFCBuild

        self.conf  = conf.MFCConf(self)
        self.user  = user.MFCUser()
        self.setup_directories()
        self.lock  = lock.MFCLock(self)
        self.args  = args.parse(self)
        self.clean = clean.MFCClean(self)
        self.build = MFCBuild(self)
        self.test  = test.MFCTest(self)
        self.run   = run.MFCRun(self)

        rich.print(common.MFC_HEADER)

        self.handle_changed_mode()

        rich.print(f"[yellow]You are currently in the [bold green]{self.lock.mode}[/bold green] mode.[/yellow]")

        if self.args["command"] == "test":
            rich.print("[bold][u]Test:[/u][/bold]")
            self.test.test()
        elif self.args["command"] == "run":
            self.run.run()
        elif self.args["command"] == "clean":
            rich.print("[bold][u]Clean:[/u][/bold]")
            self.clean.run()

        if self.args["command"] == "build":
            for target_name in [ x.name for x in self.conf.targets ]:
                if target_name in self.args["targets"]:
                    rich.print("[bold][u]Build:[/u][/bold]")
                    self.build.build_target(target_name)

        self.lock.save()
        self.user.save()

    def handle_changed_mode(self):
        def onNewMode():
            if onNewMode.flag:
                return
            
            onNewMode.flag = True

            rich.print(f'[yellow]Switching to mode [bold green]{self.args["mode"]}[/bold green]. Purging references to other modes...[/yellow]')
            
            # Update mode in mfc.user.yaml
            self.lock.mode = self.args["mode"]
            self.lock.save()

            for mode in self.user.modes:
                mode: user.Mode

                if mode.name == self.lock.mode:
                    return

                common.delete_directory_recursive(self.build.get_mode_base_path(mode.name))
        
        onNewMode.flag = False

        for idx, entry in enumerate(self.lock.targets):
            entry: lock.LockTargetHolder

            if entry.target.common_mode == None and entry.metadata.mode != self.args["mode"]:
                onNewMode()
                
                del self.lock.targets[idx]
            
            self.lock.save()

    def setup_directories(self):
        common.create_directory(common.MFC_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for mode in [ mode.name for mode in self.user.modes ] + ["common"]:
                common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        common.create_directory(f"{common.MFC_SUBDIR}/{mode}/{d}/{build_subdir}")

def main():
    mfc = MFCState()

if __name__ == "__main__":
    try:
        main()
    except common.MFCException as exc:
#        traceback.print_exc()
        rich.print(f"[red]> {str(exc)}[/red]")
        exit(1)
    except KeyboardInterrupt as exc:
        exit(1)
