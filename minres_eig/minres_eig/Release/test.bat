@echo off
cls
set /a size=40

:loop
set /a l=1

:iloop
echo [%size%] loop:%l%

echo %size% 10>val
minres_eig.exe<val>%size%_%l%.txt

if "%l%"=="5" goto next

set /a l=l+1
goto iloop

:next
if "%size%"=="120" goto end

set /a size=size+20
goto loop

:end
pause