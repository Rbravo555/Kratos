# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix

# Include any dependencies generated for this target.
include example/CMakeFiles/reentrant.dir/depend.make

# Include the progress variables for this target.
include example/CMakeFiles/reentrant.dir/progress.make

# Include the compile flags for this target's objects.
include example/CMakeFiles/reentrant.dir/flags.make

example/CMakeFiles/reentrant.dir/reentrant.c.o: example/CMakeFiles/reentrant.dir/flags.make
example/CMakeFiles/reentrant.dir/reentrant.c.o: example/reentrant.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object example/CMakeFiles/reentrant.dir/reentrant.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/reentrant.dir/reentrant.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/reentrant.c

example/CMakeFiles/reentrant.dir/reentrant.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/reentrant.dir/reentrant.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/reentrant.c > CMakeFiles/reentrant.dir/reentrant.c.i

example/CMakeFiles/reentrant.dir/reentrant.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/reentrant.dir/reentrant.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/reentrant.c -o CMakeFiles/reentrant.dir/reentrant.c.s

example/CMakeFiles/reentrant.dir/reentrant.c.o.requires:

.PHONY : example/CMakeFiles/reentrant.dir/reentrant.c.o.requires

example/CMakeFiles/reentrant.dir/reentrant.c.o.provides: example/CMakeFiles/reentrant.dir/reentrant.c.o.requires
	$(MAKE) -f example/CMakeFiles/reentrant.dir/build.make example/CMakeFiles/reentrant.dir/reentrant.c.o.provides.build
.PHONY : example/CMakeFiles/reentrant.dir/reentrant.c.o.provides

example/CMakeFiles/reentrant.dir/reentrant.c.o.provides.build: example/CMakeFiles/reentrant.dir/reentrant.c.o


# Object files for target reentrant
reentrant_OBJECTS = \
"CMakeFiles/reentrant.dir/reentrant.c.o"

# External object files for target reentrant
reentrant_EXTERNAL_OBJECTS =

example/reentrant: example/CMakeFiles/reentrant.dir/reentrant.c.o
example/reentrant: example/CMakeFiles/reentrant.dir/build.make
example/reentrant: /usr/lib/x86_64-linux-gnu/libhwloc.so
example/reentrant: /home/josep_maria/Kratos/libs/libscotch.a
example/reentrant: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
example/reentrant: /usr/lib/x86_64-linux-gnu/libpthread.so
example/reentrant: /usr/lib/x86_64-linux-gnu/libz.so
example/reentrant: /usr/lib/x86_64-linux-gnu/libm.so
example/reentrant: /usr/lib/x86_64-linux-gnu/librt.so
example/reentrant: libpastix.a
example/reentrant: /usr/lib/x86_64-linux-gnu/libopenblas.so
example/reentrant: spm/libspm.a
example/reentrant: kernels/libpastix_kernels.a
example/reentrant: /usr/lib/x86_64-linux-gnu/libhwloc.so
example/reentrant: /home/josep_maria/Kratos/libs/libscotch.a
example/reentrant: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
example/reentrant: /usr/lib/x86_64-linux-gnu/libpthread.so
example/reentrant: /usr/lib/x86_64-linux-gnu/libz.so
example/reentrant: /usr/lib/x86_64-linux-gnu/liblapacke.so
example/reentrant: /usr/lib/x86_64-linux-gnu/libm.so
example/reentrant: /usr/lib/x86_64-linux-gnu/librt.so
example/reentrant: /usr/lib/x86_64-linux-gnu/libopenblas.so
example/reentrant: example/CMakeFiles/reentrant.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable reentrant"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/reentrant.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
example/CMakeFiles/reentrant.dir/build: example/reentrant

.PHONY : example/CMakeFiles/reentrant.dir/build

example/CMakeFiles/reentrant.dir/requires: example/CMakeFiles/reentrant.dir/reentrant.c.o.requires

.PHONY : example/CMakeFiles/reentrant.dir/requires

example/CMakeFiles/reentrant.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && $(CMAKE_COMMAND) -P CMakeFiles/reentrant.dir/cmake_clean.cmake
.PHONY : example/CMakeFiles/reentrant.dir/clean

example/CMakeFiles/reentrant.dir/depend:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/CMakeFiles/reentrant.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : example/CMakeFiles/reentrant.dir/depend

