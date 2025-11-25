#!/bin/zsh
if ! ./configure.sh; then
    echo "Configure failed"
    exit 1
fi

if [[ $# -eq 0 ]]; then
    config="RelWithDebInfo"
else
    config="$1"
fi

cmake --build build --config "${config}" || {
    echo "Build failed"
    exit 1
}

echo "Build succeeded"
