@set PATH=%~dp0\w64devkit\bin;%PATH%
@set HOME=%~dp0
@busybox sh --login -i -c "cd '%~dp0' && ./w64devkit/compile.sh"