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
include example/CMakeFiles/dump_rank.dir/depend.make

# Include the progress variables for this target.
include example/CMakeFiles/dump_rank.dir/progress.make

# Include the compile flags for this target's objects.
include example/CMakeFiles/dump_rank.dir/flags.make

example/CMakeFiles/dump_rank.dir/dump_rank.c.o: example/CMakeFiles/dump_rank.dir/flags.make
example/CMakeFiles/dump_rank.dir/dump_rank.c.o: example/dump_rank.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object example/CMakeFiles/dump_rank.dir/dump_rank.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/dump_rank.dir/dump_rank.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/dump_rank.c

example/CMakeFiles/dump_rank.dir/dump_rank.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/dump_rank.dir/dump_rank.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/dump_rank.c > CMakeFiles/dump_rank.dir/dump_rank.c.i

example/CMakeFiles/dump_rank.dir/dump_rank.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/dump_rank.dir/dump_rank.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/dump_rank.c -o CMakeFiles/dump_rank.dir/dump_rank.c.s

example/CMakeFiles/dump_rank.dir/dump_rank.c.o.requires:

.PHONY : example/CMakeFiles/dump_rank.dir/dump_rank.c.o.requires

example/CMakeFiles/dump_rank.dir/dump_rank.c.o.provides: example/CMakeFiles/dump_rank.dir/dump_rank.c.o.requires
	$(MAKE) -f example/CMakeFiles/dump_rank.dir/build.make example/CMakeFiles/dump_rank.dir/dump_rank.c.o.provides.build
.PHONY : example/CMakeFiles/dump_rank.dir/dump_rank.c.o.provides

example/CMakeFiles/dump_rank.dir/dump_rank.c.o.provides.build: example/CMakeFiles/dump_rank.dir/dump_rank.c.o


# Object files for target dump_rank
dump_rank_OBJECTS = \
"CMakeFiles/dump_rank.dir/dump_rank.c.o"

# External object files for target dump_rank
dump_rank_EXTERNAL_OBJECTS =

example/dump_rank: example/CMakeFiles/dump_rank.dir/dump_rank.c.o
example/dump_rank: example/CMakeFiles/dump_rank.dir/build.make
example/dump_rank: /usr/lib/x86_64-linux-gnu/libhwloc.so
example/dump_rank: /home/josep_maria/Kratos/libs/libscotch.a
example/dump_rank: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/libpthread.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/libz.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/libm.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/librt.so
example/dump_rank: libpastix.a
example/dump_rank: /usr/lib/x86_64-linux-gnu/libopenblas.so
example/dump_rank: spm/libspm.a
example/dump_rank: kernels/libpastix_kernels.a
example/dump_rank: /usr/lib/x86_64-linux-gnu/libhwloc.so
example/dump_rank: /home/josep_maria/Kratos/libs/libscotch.a
example/dump_rank: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/libpthread.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/libz.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/liblapacke.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/libm.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/librt.so
example/dump_rank: /usr/lib/x86_64-linux-gnu/libopenblas.so
example/dump_rank: example/CMakeFiles/dump_rank.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable dump_rank"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dump_rank.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
example/CMakeFiles/dump_rank.dir/build: example/dump_rank

.PHONY : example/CMakeFiles/dump_rank.dir/build

example/CMakeFiles/dump_rank.dir/requires: example/CMakeFiles/dump_rank.dir/dump_rank.c.o.requires

.PHONY : example/CMakeFiles/dump_rank.dir/requires

example/CMakeFiles/dump_rank.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && $(CMAKE_COMMAND) -P CMakeFiles/dump_rank.dir/cmake_clean.cmake
.PHONY : example/CMakeFiles/dump_rank.dir/clean

example/CMakeFiles/dump_rank.dir/depend:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/CMakeFiles/dump_rank.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : example/CMakeFiles/dump_rank.dir/depend

