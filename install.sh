#!/bin/zsh
if [[ $# -eq 0 ]]; then
    config="RelWithDebInfo"
else
    config="$1"
fi

if ! ./build.sh ${config}; then
    echo "Build failed"
    exit 1
fi

rm -fr qATE.app
cmake --install build --config ${config} || {
    echo "Install failed"
    exit 1
}

echo "Install succeeded"
