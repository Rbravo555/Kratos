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
include test/CMakeFiles/c_ge2lr_tests.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/c_ge2lr_tests.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/c_ge2lr_tests.dir/flags.make

test/c_ge2lr_tests.c: test/z_ge2lr_tests.c
test/c_ge2lr_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/c_ge2lr_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating c_ge2lr_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f c_ge2lr_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_ge2lr_tests.c -p c -P ./ && chmod a-w c_ge2lr_tests.c

test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o: test/CMakeFiles/c_ge2lr_tests.dir/flags.make
test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o: test/c_ge2lr_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -o CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_ge2lr_tests.c

test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_ge2lr_tests.c > CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.i

test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_ge2lr_tests.c -o CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.s

test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.requires:

.PHONY : test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.requires

test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.provides: test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/c_ge2lr_tests.dir/build.make test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.provides.build
.PHONY : test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.provides

test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.provides.build: test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o


# Object files for target c_ge2lr_tests
c_ge2lr_tests_OBJECTS = \
"CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o"

# External object files for target c_ge2lr_tests
c_ge2lr_tests_EXTERNAL_OBJECTS =

test/c_ge2lr_tests: test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o
test/c_ge2lr_tests: test/CMakeFiles/c_ge2lr_tests.dir/build.make
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libhwloc.so
test/c_ge2lr_tests: /home/josep_maria/Kratos/libs/libscotch.a
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libpthread.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libz.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/librt.so
test/c_ge2lr_tests: test/libpastix_tests.a
test/c_ge2lr_tests: kernels/libpastix_kernels.a
test/c_ge2lr_tests: libpastix.a
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/liblapacke.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libtmglib.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/c_ge2lr_tests: kernels/libpastix_kernels.a
test/c_ge2lr_tests: spm/libspm.a
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libhwloc.so
test/c_ge2lr_tests: /home/josep_maria/Kratos/libs/libscotch.a
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libscotcherrexit.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libpthread.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libz.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/librt.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/liblapacke.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libm.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/librt.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/liblapacke.so
test/c_ge2lr_tests: /usr/lib/x86_64-linux-gnu/libopenblas.so
test/c_ge2lr_tests: test/CMakeFiles/c_ge2lr_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable c_ge2lr_tests"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/c_ge2lr_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/c_ge2lr_tests.dir/build: test/c_ge2lr_tests

.PHONY : test/CMakeFiles/c_ge2lr_tests.dir/build

test/CMakeFiles/c_ge2lr_tests.dir/requires: test/CMakeFiles/c_ge2lr_tests.dir/c_ge2lr_tests.c.o.requires

.PHONY : test/CMakeFiles/c_ge2lr_tests.dir/requires

test/CMakeFiles/c_ge2lr_tests.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && $(CMAKE_COMMAND) -P CMakeFiles/c_ge2lr_tests.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/c_ge2lr_tests.dir/clean

test/CMakeFiles/c_ge2lr_tests.dir/depend: test/c_ge2lr_tests.c
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/CMakeFiles/c_ge2lr_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/c_ge2lr_tests.dir/depend

