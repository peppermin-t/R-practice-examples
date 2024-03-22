
@echo off
for /l %%I in (1, 1, 100) do (
	echo "trial%%I:"
	Rscript.exe G03.R %%I
)
pause

