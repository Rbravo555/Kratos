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
include spm/tests/CMakeFiles/spm_dof_norm_tests.dir/depend.make

# Include the progress variables for this target.
include spm/tests/CMakeFiles/spm_dof_norm_tests.dir/progress.make

# Include the compile flags for this target's objects.
include spm/tests/CMakeFiles/spm_dof_norm_tests.dir/flags.make

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o: spm/tests/CMakeFiles/spm_dof_norm_tests.dir/flags.make
spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o: spm/tests/spm_dof_norm_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests/spm_dof_norm_tests.c

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests/spm_dof_norm_tests.c > CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.i

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests/spm_dof_norm_tests.c -o CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.s

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.requires:

.PHONY : spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.requires

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.provides: spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.requires
	$(MAKE) -f spm/tests/CMakeFiles/spm_dof_norm_tests.dir/build.make spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.provides.build
.PHONY : spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.provides

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.provides.build: spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o


# Object files for target spm_dof_norm_tests
spm_dof_norm_tests_OBJECTS = \
"CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o"

# External object files for target spm_dof_norm_tests
spm_dof_norm_tests_EXTERNAL_OBJECTS =

spm/tests/spm_dof_norm_tests: spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o
spm/tests/spm_dof_norm_tests: spm/tests/CMakeFiles/spm_dof_norm_tests.dir/build.make
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libhwloc.so
spm/tests/spm_dof_norm_tests: /home/josep_maria/Kratos/libs/libscotch.a
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libpthread.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libz.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libm.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/librt.so
spm/tests/spm_dof_norm_tests: spm/libspm.a
spm/tests/spm_dof_norm_tests: spm/tests/libspm_test.a
spm/tests/spm_dof_norm_tests: spm/libspm.a
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libhwloc.so
spm/tests/spm_dof_norm_tests: /home/josep_maria/Kratos/libs/libscotch.a
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libpthread.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libz.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/liblapacke.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/libm.so
spm/tests/spm_dof_norm_tests: /usr/lib/x86_64-linux-gnu/librt.so
spm/tests/spm_dof_norm_tests: spm/tests/CMakeFiles/spm_dof_norm_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable spm_dof_norm_tests"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/spm_dof_norm_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
spm/tests/CMakeFiles/spm_dof_norm_tests.dir/build: spm/tests/spm_dof_norm_tests

.PHONY : spm/tests/CMakeFiles/spm_dof_norm_tests.dir/build

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/requires: spm/tests/CMakeFiles/spm_dof_norm_tests.dir/spm_dof_norm_tests.c.o.requires

.PHONY : spm/tests/CMakeFiles/spm_dof_norm_tests.dir/requires

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests && $(CMAKE_COMMAND) -P CMakeFiles/spm_dof_norm_tests.dir/cmake_clean.cmake
.PHONY : spm/tests/CMakeFiles/spm_dof_norm_tests.dir/clean

spm/tests/CMakeFiles/spm_dof_norm_tests.dir/depend:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/tests/CMakeFiles/spm_dof_norm_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : spm/tests/CMakeFiles/spm_dof_norm_tests.dir/depend

