#!/bin/zsh
local script_dir=$(dirname $(readlink -f "$0"))
source ${script_dir}/qt_env_macos.sh
cmake $@
