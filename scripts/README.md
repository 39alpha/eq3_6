# Scripts

The shell script eq36cfg configures the user environment to run this version of EQ3/6 by
making necessary path changes and setting required enviroment variables. It can be run from
a user's .cshrc file. This shell script is functionally equivalent to the eq36cfg.bat batch
file used on Windows systems.

The shell scripts runeq36, runeqpt, and xcif36 respective run EQ3NR and EQ6, RUNEQPT, and
XCON3 and XCON6. To use runeq36 to run EQ3NR, use a symbolic link named runeq3. To use
it to run EQ6, use a symbolic link named runeq6. Similarly, use symbolic links to xcif36
named xcif3 and xcif6 to run, respectively, XCON3 and XCON6. The equivalent Fortran programs
found in the icodes folder may be used instead. They are the programs used on Windows
systems.

These scripts have been tested on a Linux system (aztec.llnl.gov) and found to work
satisfactorily. They were originally written for use on Solaris systems and are thought
likely to work on modern Solaris systems.

These scripts are all C-shell type.
