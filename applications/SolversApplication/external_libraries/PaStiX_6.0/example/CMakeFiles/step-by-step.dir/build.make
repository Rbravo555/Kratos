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
include example/CMakeFiles/step-by-step.dir/depend.make

# Include the progress variables for this target.
include example/CMakeFiles/step-by-step.dir/progress.make

# Include the compile flags for this target's objects.
include example/CMakeFiles/step-by-step.dir/flags.make

example/CMakeFiles/step-by-step.dir/step-by-step.c.o: example/CMakeFiles/step-by-step.dir/flags.make
example/CMakeFiles/step-by-step.dir/step-by-step.c.o: example/step-by-step.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object example/CMakeFiles/step-by-step.dir/step-by-step.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/step-by-step.dir/step-by-step.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/step-by-step.c

example/CMakeFiles/step-by-step.dir/step-by-step.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/step-by-step.dir/step-by-step.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/step-by-step.c > CMakeFiles/step-by-step.dir/step-by-step.c.i

example/CMakeFiles/step-by-step.dir/step-by-step.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/step-by-step.dir/step-by-step.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/step-by-step.c -o CMakeFiles/step-by-step.dir/step-by-step.c.s

example/CMakeFiles/step-by-step.dir/step-by-step.c.o.requires:

.PHONY : example/CMakeFiles/step-by-step.dir/step-by-step.c.o.requires

example/CMakeFiles/step-by-step.dir/step-by-step.c.o.provides: example/CMakeFiles/step-by-step.dir/step-by-step.c.o.requires
	$(MAKE) -f example/CMakeFiles/step-by-step.dir/build.make example/CMakeFiles/step-by-step.dir/step-by-step.c.o.provides.build
.PHONY : example/CMakeFiles/step-by-step.dir/step-by-step.c.o.provides

example/CMakeFiles/step-by-step.dir/step-by-step.c.o.provides.build: example/CMakeFiles/step-by-step.dir/step-by-step.c.o


# Object files for target step-by-step
step__by__step_OBJECTS = \
"CMakeFiles/step-by-step.dir/step-by-step.c.o"

# External object files for target step-by-step
step__by__step_EXTERNAL_OBJECTS =

example/step-by-step: example/CMakeFiles/step-by-step.dir/step-by-step.c.o
example/step-by-step: example/CMakeFiles/step-by-step.dir/build.make
example/step-by-step: /usr/lib/x86_64-linux-gnu/libhwloc.so
example/step-by-step: /home/josep_maria/Kratos/libs/libscotch.a
example/step-by-step: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/libpthread.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/libz.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/libm.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/librt.so
example/step-by-step: libpastix.a
example/step-by-step: /usr/lib/x86_64-linux-gnu/libopenblas.so
example/step-by-step: spm/libspm.a
example/step-by-step: kernels/libpastix_kernels.a
example/step-by-step: /usr/lib/x86_64-linux-gnu/libhwloc.so
example/step-by-step: /home/josep_maria/Kratos/libs/libscotch.a
example/step-by-step: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/libpthread.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/libz.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/liblapacke.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/libm.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/librt.so
example/step-by-step: /usr/lib/x86_64-linux-gnu/libopenblas.so
example/step-by-step: example/CMakeFiles/step-by-step.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable step-by-step"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/step-by-step.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
example/CMakeFiles/step-by-step.dir/build: example/step-by-step

.PHONY : example/CMakeFiles/step-by-step.dir/build

example/CMakeFiles/step-by-step.dir/requires: example/CMakeFiles/step-by-step.dir/step-by-step.c.o.requires

.PHONY : example/CMakeFiles/step-by-step.dir/requires

example/CMakeFiles/step-by-step.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example && $(CMAKE_COMMAND) -P CMakeFiles/step-by-step.dir/cmake_clean.cmake
.PHONY : example/CMakeFiles/step-by-step.dir/clean

example/CMakeFiles/step-by-step.dir/depend:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/example/CMakeFiles/step-by-step.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : example/CMakeFiles/step-by-step.dir/depend

