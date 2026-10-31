nwdynam: C source of the native plugin
======================================

nwdynam runs its heavy computation in a compiled Stata plugin (dynam_sim). The package
installs the compiled plugins for macOS, Linux and Windows; this file comes with
the C source they were built from. nwdynam is free software under the GNU General
Public License, version 3 (file LICENSE in the repository).

Where to get the source
  * In Stata:  net get nwdynam, from(https://raw.githubusercontent.com/thomas-u-grund/nwdynam/main)
    copies the files into the current directory. net get puts all files of a
    package into one folder; the build needs the layout of the repository (below).
  * On GitHub: https://github.com/thomas-u-grund/nwdynam, folder native/.

Files
  dynam_sim.c
  stplugin.c, stplugin.h   StataCorp's Stata Plugin Interface 3.0 (https://www.stata.com/plugins/),
                           redistributed unchanged; not part of nwdynam's own code
  makefile                 the build rules (Makefile in the repository; net get writes lower case names)

Layout needed for the build (as in the repository)
  native/dynam_sim.c native/Makefile
  native/spi/stplugin.c  native/spi/stplugin.h

Build
  macOS (universal binary):  cd native && make macos
  Linux:                     cd native && make unix
  Windows (Visual Studio x64 prompt, in native/):
      cl /LD /DSYSTEM=STWIN32 /I spi dynam_sim.c spi\stplugin.c /Fe:..\lib\plugins\windows\dynam_sim_windows.plugin
  The binaries are written to lib/plugins/<platform>/ next to the .ado files.
