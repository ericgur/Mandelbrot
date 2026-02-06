Get-ChildItem -Path src -Include *.cpp,*.hpp,*.h -File -Recurse | ForEach-Object { clang-format -i $_.FullName }
