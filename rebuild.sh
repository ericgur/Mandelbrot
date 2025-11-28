#!/bin/zsh
rm -fr build
./configure.sh || {
    echo "Configure failed!"
    exit 1
}

./build.sh "$@"
