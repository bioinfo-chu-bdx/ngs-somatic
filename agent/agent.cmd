@echo off
set CURRENT_WOR_DIR=%cd%
setlocal enabledelayedexpansion
set CMD_LINE_ARGS=%*

if not "%ATOOLS_HOME%" == "" goto OkHome
set ATOOLS_HOME=%~dp0%

:OkHome
cd /d %ATOOLS_HOME%

if not "%JAVA_HOME%" == "" if exist "%JAVA_HOME%\bin\java.exe" goto OkJavaHome
set FOUND=
for %%e in (%PATHEXT%) do (
  for %%X in (java%%e) do (
    if not defined FOUND (
      set FOUND=%%~$PATH:X
    )
  )
)
if defined FOUND goto OkJava
echo "No JAVA_HOME defined and no Java found in PATH. Please intstall Java 8 or higher to use AGeNT"
exit 0


:OkJava
set _JAVA="java"


:OkJavaHome
set _JAVA="%JAVA_HOME%\bin\java.exe"

for /F "tokens=1" %%a in ('dir "lib\Agent*.jar" /b') do set trimmerJarName=%%a

%_JAVA% -jar "lib/%trimmerJarName%" "%CURRENT_WOR_DIR%" -home:%ATOOLS_HOME% %*