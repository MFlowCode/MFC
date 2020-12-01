# How to get MFC up and running on Richardson

* Log in
    * $ ssh username@richardson.caltech.edu
    * Make this an alias! e.g. $ richardson
        * Append 
        alias richardson='ssh -Y <your_userid>@richardson.caltech.edu'
        RICHARDSON="<user name>@richardson.caltech.edu"
        to your ~/.bash_profile  
        on your local computer
    * Add your private key to your Richardson keychain
        * $ cat ~/.ssh/id_rsa.pub | ssh <your_userid>@richardson.caltech.edu 'cat >> ~/.ssh/authorized_keys'

=======
        * [Generate](https://docs.joyent.com/public-cloud/getting-started/ssh-keys/generating-an-ssh-key-manually/manually-generating-your-ssh-key-in-mac-os-x) your private key if you don't have one 
        * $ cat ~/.ssh/id_rsa.pub | ssh <your_userid>@richardson.caltech.edu 'cat >> ~/.ssh/authorized_keys'

* Silo/HDF5, Visit, Paraview on Richardson
    * Follow the `visit_paraview_richardson.pdf` document  
    * Step 15. under 'visit on richardson' should use the following URL:  http://visit.ilight.com/svn/visit/branches/2.9RC/src/svn_bin/bv_support/

>>>>>>> qbmm
* Clone MFC, e.g.
    * $ git clone https://github.com/ComputationalFlowPhysics/MFC_private.git 
    * $ git status
        * shows what branch you are on
    * $ git checkout <branch>
        * switch to branch <branch> 

* Compile MFC
    * Get open-mpi
        * Show currently loaded modules
            * $ module list  
        * Show possible module files
            * $ module avail
        * Load open-mpi
            * For example: $ module load openmpi-1.8/gcc 
        * Add this to your .bashrc so it happens at login time
            * Append 
            module load openmpi-1.8/gcc
            to your ~/.bashrc
            on Richardson
    * Install FFTW
        * $ cd installers
        * $ ./install_fftw.sh
    * Install LAPACK [if QBMM]
        * $ cd installers
        * $ ./install_lapack.sh
    * Make MFC
        * $ cd ../    [or to MFC directory]
        * make
        * make test

* Moving files from local computer to Richardson [While on local computer. Can't do other way around]
    * $ rsync -rvz <local files> <user name>@richardson.caltech.edu:<richardson location>
    * or $ rsync -rvz <local files> $RICHARDSON:<richardson location>
    * Or from Richardson to local computer
        * $ rsync -rvz <user name>@richardson.caltech.edu:<richardson location> <local files>
        * or $ rsync -rvz $RICHARDSON:<richardson location> <local files>
    * Can also use GUI-based tools like FileZilla 

* How to visualize using gnuplot on Richardson
    * You need XQuartz on your MacOS [or X11 on Linux]
        * Install at https://dl.bintray.com/xquartz/downloads/XQuartz-2.7.8.dmg
        * Make sure you logged into Richardson using ssh -Y [or ssh -X] 
        * Check if your X server is working by logging into Richardson and issuing
            $ xclock
            * I am aware this version is deprecated, but it appears that 2.7.11 does not work in some cases
        * Make sure you logged into Richardson using ssh -Y [or ssh -X] 
        * Check if your X server is working by logging into Richardson and issuing
            $ xclock
        * IF YOU USE MacOS CATALINA!
            * Permissions can be an issue
            * Go to System preferences -> Security and Privacy -> Privacy -> Full Disk Access 
            * Add your terminal (e.g. Terminal, iTerm2, etc.) via the + button
            * Add XQuartz (launchd_startx) and the app (XQuartz)

    * Copy contents of https://raw.githubusercontent.com/sbryngelson/dotfiles/master/.gnuplotrc_x11
    into your ~/.bashrc
        * now you can use 
        $ myplot cons.1.0000000
        for, e.g., conservative and primitive serial data files cons.* prim.*
        $ myplots X probe1.dat
        for probe data [variable numeber X]
        $ animate1d X cons.1.*
        to animate cons.1.* X times



