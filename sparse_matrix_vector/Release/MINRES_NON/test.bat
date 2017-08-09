@echo off
cd %~dp0
cls
set /a size=50

:loop
set /a l=1

:iloop
echo [%size%] loop:%l%

echo %size% 05> val
sparse_matrix_vector.exe<val

mkdir result\%l%
copy result\*.csv result\%l%\*.csv
del result\*.csv

if "%l%"=="5" goto next

set /a l=l+1
goto iloop

:next
if "%size%"=="250" goto end

set /a size=size+50
goto loop

:end