#!/bin/bash
set -e

echo "Installing required dependencies..."

# Detect OS and install packages
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Detected macOS"
    brew update
    brew install cmake protobuf boost python@3.10
elif [[ -f /etc/debian_version ]]; then
    echo "Detected Debian/Ubuntu"
    sudo apt-get update
    sudo apt-get install -y build-essential cmake protobuf-compiler libboost-all-dev python3 python3-pip
elif [[ -f /etc/redhat-release ]]; then
    echo "Detected RedHat/CentOS"
    sudo yum install -y epel-release
    sudo yum install -y cmake3 protobuf-compiler boost-devel python3
else
    echo "Unsupported OS: $OSTYPE"
    exit 1
fi

echo "Dependencies installed successfully."
