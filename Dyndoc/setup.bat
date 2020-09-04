@echo off
SET PATHROOT=C:\Program Files\R\
echo Locating path of R...
echo.
if not exist "%PATHROOT%" goto:NO_R
for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do (
echo Found %%r
echo shell "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do
echo All set!
goto:DONE
)
:NO_R
echo R is not installed in your system.
echo.
echo Download it from https://cran.r-project.org/bin/windows/base/
echo Install it and re-run this script
:DONE
echo.
pause