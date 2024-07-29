#!/bin/bash

if ! command -v docker > /dev/null 2>&1; then
    error "$MAGENTA""Docker$COLOR_RESET is not installed."

    exit 1
fi

log "Running in$MAGENTA Docker$COLOR_RESET mode."

if [ -t 1 ]; then
    dockerintopts='--interactive --tty'
fi

__docker_run() {
    docker run $dockerintopts --rm --workdir /home/me/MFC            \
               --mount type=bind,source="$(pwd)",target=/home/me/MFC \
               sbryngelson/mfc:latest $@
}

__docker_run sudo chown -R me:me /home/me/MFC
if (($?)); then
    error "Docker: Failed to set directory permissions on MFC mount.."

    exit 1
fi

__docker_run $@
if (($?)); then
    error "Error running Docker container with $@."

    exit 1
fi
