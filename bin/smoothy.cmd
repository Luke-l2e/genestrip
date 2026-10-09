@echo off
setlocal

rem Starts the smoothy module natively on Windows, like bin/smoothy.sh does under Linux/WSL.
rem It runs lib\smoothy.jar, which is built by `mvn -pl smoothy -am package`.
rem Heap as in bin/smoothy.sh; override per run with GS_XMX if a machine has more or less.
if not defined GS_XMX set GS_XMX=56g

java -Xmx%GS_XMX% -jar "%~dp0..\lib\smoothy.jar" -d "%~dp0..\data" %*
