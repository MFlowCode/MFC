#!/usr/bin/env bash

if [ ! -f ./bootstrap/mfc.py ]; then
    echo "[mfc.sh] Error: You must call this script from within MFC's root folder."
    exit 1
fi

which python3 > /dev/null 2>&1
if (($?)); then
    echo "[mfc.sh] Error: Couldn't find Python3. Please ensure it is discoverable."
    exit 1
fi

python3 -c 'print("")' > /dev/null
if (($?)); then
    echo "[mfc.sh] Error: Python3 is present but can't execute a simple program. Please ensure that python3 is working."
    exit 1
fi

chmod +x ./bootstrap/mfc.py
if (($?)); then
    echo "[mfc.sh] Error: Failed to set ./bootstrap/mfc.py as executable."
    exit 1
fi

cd bootstrap && ./mfc.py $@
exit $?
