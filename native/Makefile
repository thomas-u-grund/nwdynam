# Native (C) Stata plugins of the nwdynam package: dynam_sim.
# Output: ../lib/plugins/<platform>/<plugin>_<platform>.plugin (the platform suffix matters for `net install`,
# which flattens all files into one folder). spi/ holds StataCorp's plugin interface (SPI 3.0, redistributed as
# published at https://www.stata.com/plugins/); it is not this package's own code.
#   make macos    universal (x86_64 + arm64) binaries, run on macOS
#   make unix     Linux (run on a Linux host or the CI runner)
# Windows (Visual Studio x64 prompt): cl /LD /DSYSTEM=STWIN32 /I spi <plugin>.c spi\stplugin.c /Fe:..\lib\plugins\windows\<plugin>_windows.plugin

CC = clang
CFLAGS = -O2 -Wall -Ispi

.PHONY: macos unix clean macos-dynam_sim unix-dynam_sim

macos: macos-dynam_sim

unix: unix-dynam_sim

macos-dynam_sim:
	mkdir -p ../lib/plugins/macos
	$(CC) -bundle -o dynam_sim.x86_64 dynam_sim.c spi/stplugin.c \
		-DSYSTEM=APPLEMAC -target x86_64-apple-macos10.11 $(CFLAGS)
	$(CC) -bundle -o dynam_sim.arm64 dynam_sim.c spi/stplugin.c \
		-DSYSTEM=APPLEMAC -target arm64-apple-macos11 $(CFLAGS)
	lipo -create -output ../lib/plugins/macos/dynam_sim_macos.plugin dynam_sim.x86_64 dynam_sim.arm64
	codesign -s - ../lib/plugins/macos/dynam_sim_macos.plugin
	rm -f dynam_sim.x86_64 dynam_sim.arm64
	@echo "Built ../lib/plugins/macos/dynam_sim_macos.plugin"

unix-dynam_sim:
	mkdir -p ../lib/plugins/unix
	gcc -shared -fPIC -DSYSTEM=OPUNIX $(CFLAGS) -o ../lib/plugins/unix/dynam_sim_unix.plugin \
		dynam_sim.c spi/stplugin.c -ldl
	@echo "Built ../lib/plugins/unix/dynam_sim_unix.plugin"

clean:
	rm -f dynam_sim.x86_64 dynam_sim.arm64 ../lib/plugins/macos/dynam_sim_macos.plugin ../lib/plugins/unix/dynam_sim_unix.plugin
