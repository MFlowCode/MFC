* New branches cannot be made on MFlowCode/MFC, they are made on forks
* PRs:
    * made using AI tools like Claude Code and Codex should say so. 
    * are made from those MFC forks
    * that change CFD result need verification PR is correct
    * follow template
    * that break a feature but promise a followup PR to fix it are rejected
* Commands:
    * MFC should almost always build and run using the ./mfc.sh command
    * Running mfc.sh commands can create a lock file in build/ that is sticky, be careful
* Programming and Design:
    * New code should follow the DRY principle and also make side-effect code DRY as well
    * Comments should be as short as possible without sacrificing value
    * GPU macros should follow existing the source's GPU macro principles and patterns
    * Functions/subroutines/modules shorter is better while being correctness, fast, and separating concerns 
