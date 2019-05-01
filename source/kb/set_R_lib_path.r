# How to set shared lib folder for all R versions
# April 26, 2019

# Windows

1. Create a new folder, eg. C:\users\abc\Rlib. Optionally, copy old packages to this folder.

2. open any text editor and save "R_LIBS_USER=C:\users\abc\Rlib" in a file
named .Renviron under your HOME directory (usually C:\users\abc\Documents)
# the correct path can be obtained by R command path.expand("~")

# or go to  Start --> Control Panel --> User Accounts --> Change my environmental variables,
# add a new entry, R_LIBS_USER - C:\users\abc\Rlib ( note backslashes in Windows)

3. Within R, run
update.packages(checkBuilt = TRUE, ask = FALSE)



# Linux

1. Create a new folder, eg. /home/abc/R/library and optionally copy old packages to this folder
2. in terminal run:
  echo "R_LIBS_USER=/home/abc/R/library" >/home/abc/.Renviron

3. Within R, run
update.packages(checkBuilt = TRUE, ask = FALSE)
