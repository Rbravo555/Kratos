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
include test/CMakeFiles/d_lrmm_tests.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/d_lrmm_tests.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/d_lrmm_tests.dir/flags.make

test/d_lrmm_tests.c: test/z_lrmm_tests.c
test/d_lrmm_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/d_lrmm_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating d_lrmm_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f d_lrmm_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_lrmm_tests.c -p d -P ./ && chmod a-w d_lrmm_tests.c

test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o: test/CMakeFiles/d_lrmm_tests.dir/flags.make
test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o: test/d_lrmm_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -o CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_lrmm_tests.c

test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_lrmm_tests.c > CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.i

test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_lrmm_tests.c -o CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.s

test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.requires:

.PHONY : test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.requires

test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.provides: test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/d_lrmm_tests.dir/build.make test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.provides.build
.PHONY : test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.provides

test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.provides.build: test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o


# Object files for target d_lrmm_tests
d_lrmm_tests_OBJECTS = \
"CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o"

# External object files for target d_lrmm_tests
d_lrmm_tests_EXTERNAL_OBJECTS =

test/d_lrmm_tests: test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o
test/d_lrmm_tests: test/CMakeFiles/d_lrmm_tests.dir/build.make
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libhwloc.so
test/d_lrmm_tests: /home/josep_maria/Kratos/libs/libscotch.a
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libpthread.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libz.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/librt.so
test/d_lrmm_tests: test/libpastix_tests.a
test/d_lrmm_tests: kernels/libpastix_kernels.a
test/d_lrmm_tests: libpastix.a
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/liblapacke.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libtmglib.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/d_lrmm_tests: kernels/libpastix_kernels.a
test/d_lrmm_tests: spm/libspm.a
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libhwloc.so
test/d_lrmm_tests: /home/josep_maria/Kratos/libs/libscotch.a
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libpthread.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libz.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/librt.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/liblapacke.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/librt.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/liblapacke.so
test/d_lrmm_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/d_lrmm_tests: test/CMakeFiles/d_lrmm_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable d_lrmm_tests"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/d_lrmm_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/d_lrmm_tests.dir/build: test/d_lrmm_tests

.PHONY : test/CMakeFiles/d_lrmm_tests.dir/build

test/CMakeFiles/d_lrmm_tests.dir/requires: test/CMakeFiles/d_lrmm_tests.dir/d_lrmm_tests.c.o.requires

.PHONY : test/CMakeFiles/d_lrmm_tests.dir/requires

test/CMakeFiles/d_lrmm_tests.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && $(CMAKE_COMMAND) -P CMakeFiles/d_lrmm_tests.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/d_lrmm_tests.dir/clean

test/CMakeFiles/d_lrmm_tests.dir/depend: test/d_lrmm_tests.c
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/CMakeFiles/d_lrmm_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/d_lrmm_tests.dir/depend

