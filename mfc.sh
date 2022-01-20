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

#python3 -m venv ./venv
#if (($?)); then
#    echo "[mfc.sh] Error: Failed to create a Python virtual environment."
#    exit 1
#fi

wget -O ./bootstrap/get-pip.py https://bootstrap.pypa.io/get-pip.py
if (($?)); then
    echo "[mfc.sh] Error: Couldn't download get-pip.py using wget."
    exit 1
fi

python3 ./bootstrap/get-pip.py --user
if (($?)); then
    echo "[mfc.sh] Error: Coudln't install pip with get-pip.py."
    exit 1
fi

#source ./venv/bin/activate
#if (($?)); then
#    echo "[mfc.sh] Error: Faild to activate the Python virtual environment."
#    exit 1
#fi

python3 -c 'print("")' > /dev/null
if (($?)); then
    echo "[mfc.sh] Error: Python3 is present but can't execute a simple program. Please ensure that python3 is working."
    exit 1
fi

cd bootstrap
python3 ./mfc.py $@
cd ..

exit $?

