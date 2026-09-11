@echo off
REM LexiPair-BioRange one-command runner (Windows)
python tools\generate_datasets.py
if errorlevel 1 exit /b 1
python run_all.py
