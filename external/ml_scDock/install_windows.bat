@echo off
setlocal EnableDelayedExpansion

echo ============================================================
echo  scDock Windows Installer
echo  Installs AutoDock Vina and OpenBabel via conda
echo ============================================================
echo.

:: ---------------------------------------------------------------
:: Check conda is available
:: ---------------------------------------------------------------
where conda >nul 2>&1
if errorlevel 1 (
    echo [ERROR] conda not found on PATH.
    echo.
    echo Please install Miniconda first:
    echo   https://docs.conda.io/en/latest/miniconda.html
    echo.
    echo After installing, reopen this prompt from:
    echo   Start Menu -^> Anaconda Prompt  (or Miniconda Prompt)
    echo and run this script again.
    pause
    exit /b 1
)

echo [OK] conda found.
conda --version
echo.

:: ---------------------------------------------------------------
:: Create or reuse conda environment
:: ---------------------------------------------------------------
set ENV_NAME=scdock

conda info --envs | findstr /B "%ENV_NAME%" >nul 2>&1
if errorlevel 1 (
    echo Creating conda environment "%ENV_NAME%" with Python 3.10...
    call conda create -y -n %ENV_NAME% python=3.10 -c conda-forge
    if errorlevel 1 (
        echo [ERROR] Failed to create conda environment.
        pause
        exit /b 1
    )
    echo [OK] Environment "%ENV_NAME%" created.
) else (
    echo [OK] Conda environment "%ENV_NAME%" already exists. Skipping creation.
)
echo.

:: ---------------------------------------------------------------
:: Activate environment and install packages
:: ---------------------------------------------------------------
echo Activating environment "%ENV_NAME%"...
call conda activate %ENV_NAME%
if errorlevel 1 (
    echo [ERROR] Failed to activate conda environment.
    echo Try running this script from an Anaconda Prompt.
    pause
    exit /b 1
)
echo [OK] Environment activated.
echo.

echo Installing AutoDock Vina 1.1.2 and OpenBabel 3.1.1...
call conda install -y ^
    -c conda-forge ^
    -c bioconda ^
    autodock-vina=1.1.2 ^
    openbabel=3.1.1
if errorlevel 1 (
    echo [ERROR] conda install failed.
    echo Try running: conda update conda
    echo and then re-run this script.
    pause
    exit /b 1
)
echo [OK] AutoDock Vina and OpenBabel installed.
echo.

:: ---------------------------------------------------------------
:: Verify installations
:: ---------------------------------------------------------------
echo Verifying installations...

vina --version >nul 2>&1
if errorlevel 1 (
    echo [WARN] "vina" not found on PATH after install.
    echo        You may need to pass the full path to MATLAB (see below).
) else (
    echo [OK] AutoDock Vina:
    vina --version
)
echo.

obabel --version >nul 2>&1
if errorlevel 1 (
    echo [WARN] "obabel" not found on PATH after install.
    echo        You may need to pass the full path to MATLAB (see below).
) else (
    echo [OK] OpenBabel:
    obabel --version
)
echo.

:: ---------------------------------------------------------------
:: Locate executables and print MATLAB usage
:: ---------------------------------------------------------------
for /f "delims=" %%i in ('conda run -n %ENV_NAME% where vina 2^>nul') do set VINA_PATH=%%i
for /f "delims=" %%i in ('conda run -n %ENV_NAME% where obabel 2^>nul') do set OBABEL_PATH=%%i

echo ============================================================
echo  Installation complete.
echo ============================================================
echo.
echo Add the following paths to your MATLAB sc_dock call:
echo.
if defined VINA_PATH (
    echo   'vina_exe',   '%VINA_PATH%'
) else (
    echo   'vina_exe',   'C:\Users\%USERNAME%\.conda\envs\%ENV_NAME%\Scripts\vina.exe'
)
if defined OBABEL_PATH (
    echo   'obabel_exe', '%OBABEL_PATH%'
) else (
    echo   'obabel_exe', 'C:\Users\%USERNAME%\.conda\envs\%ENV_NAME%\Scripts\obabel.exe'
)
echo.
echo Example MATLAB call:
echo.
echo   result = sc_dock(X, g, 'compounds', "fda", ...
if defined VINA_PATH (
    echo       'vina_exe',   '%VINA_PATH%', ...
) else (
    echo       'vina_exe',   'C:\Users\%USERNAME%\.conda\envs\%ENV_NAME%\Scripts\vina.exe', ...
)
if defined OBABEL_PATH (
    echo       'obabel_exe', '%OBABEL_PATH%');
) else (
    echo       'obabel_exe', 'C:\Users\%USERNAME%\.conda\envs\%ENV_NAME%\Scripts\obabel.exe');
)
echo.
pause
