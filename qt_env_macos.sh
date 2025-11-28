#!/bin/zsh
QT_HOME="${HOME}/Qt"

# Find all directories that match the Qt version pattern X.Y.Z
versions=(${(f)"$(find "$QT_HOME" -maxdepth 1 -type d -regex '.*/[0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*' -exec basename {} \;)"})

# Sort versions using version-aware sorting (-V) and get the last one (latest)
if [[ ${#versions[@]} -gt 0 ]]; then
  QT_VER=$(printf "%s\n" "${versions[@]}" | sort -V | tail -n 1)
  echo "Found Qt version: ${QT_VER}"
  export QT_DIR="${QT_HOME}/${QT_VER}"
  export OpenMP_ROOT=$(brew --prefix libomp)
  export PATH=${QT_HOME}/Tools/CMake/CMake.app/Contents/bin:${QT_HOME}/Tools/Ninja:${QT_HOME}/${QT_VER}/macos/bin:${PATH}
  #PS1=%n@%m %1~ %#
  PS1="%# | Qt ${QT_VER} | %0~ >"
else
  echo "No Qt versions found in ${QT_HOME}"
fi

