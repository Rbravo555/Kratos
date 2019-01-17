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
include test/CMakeFiles/pastix_tests.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/pastix_tests.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/pastix_tests.dir/flags.make

test/c_lowrank_tests.c: test/z_lowrank_tests.c
test/c_lowrank_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/c_lowrank_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating c_lowrank_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f c_lowrank_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_lowrank_tests.c -p c -P ./ && chmod a-w c_lowrank_tests.c

test/s_lowrank_tests.c: test/z_lowrank_tests.c
test/s_lowrank_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/s_lowrank_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating s_lowrank_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f s_lowrank_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_lowrank_tests.c -p s -P ./ && chmod a-w s_lowrank_tests.c

test/d_lowrank_tests.c: test/z_lowrank_tests.c
test/d_lowrank_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/d_lowrank_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating d_lowrank_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f d_lowrank_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_lowrank_tests.c -p d -P ./ && chmod a-w d_lowrank_tests.c

test/c_bvec_tests.c: test/z_bvec_tests.c
test/c_bvec_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/c_bvec_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Generating c_bvec_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f c_bvec_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_bvec_tests.c -p c -P ./ && chmod a-w c_bvec_tests.c

test/s_bvec_tests.c: test/z_bvec_tests.c
test/s_bvec_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/s_bvec_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Generating s_bvec_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f s_bvec_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_bvec_tests.c -p s -P ./ && chmod a-w s_bvec_tests.c

test/d_bvec_tests.c: test/z_bvec_tests.c
test/d_bvec_tests.c: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
test/d_bvec_tests.c: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Generating d_bvec_tests.c"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/cmake -E remove -f d_bvec_tests.c && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_bvec_tests.c -p d -P ./ && chmod a-w d_bvec_tests.c

test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o: test/c_lowrank_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -o CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_lowrank_tests.c

test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_lowrank_tests.c > CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.i

test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_lowrank_tests.c -o CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.s

test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o


test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o: test/s_lowrank_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_s -UPRECISION_p -UPRECISION_d -UPRECISION_c -UPRECISION_z -o CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/s_lowrank_tests.c

test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_s -UPRECISION_p -UPRECISION_d -UPRECISION_c -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/s_lowrank_tests.c > CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.i

test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_s -UPRECISION_p -UPRECISION_d -UPRECISION_c -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/s_lowrank_tests.c -o CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.s

test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o


test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o: test/z_lowrank_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_z -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_c -o CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_lowrank_tests.c

test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_z -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_c -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_lowrank_tests.c > CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.i

test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_z -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_c -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_lowrank_tests.c -o CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.s

test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o


test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o: test/d_lowrank_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -o CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_lowrank_tests.c

test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_lowrank_tests.c > CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.i

test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_lowrank_tests.c -o CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.s

test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o


test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o: test/c_bvec_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -o CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_bvec_tests.c

test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/c_bvec_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_bvec_tests.c > CMakeFiles/pastix_tests.dir/c_bvec_tests.c.i

test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/c_bvec_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_c -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/c_bvec_tests.c -o CMakeFiles/pastix_tests.dir/c_bvec_tests.c.s

test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o


test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o: test/s_bvec_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_s -UPRECISION_p -UPRECISION_d -UPRECISION_c -UPRECISION_z -o CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/s_bvec_tests.c

test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/s_bvec_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_s -UPRECISION_p -UPRECISION_d -UPRECISION_c -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/s_bvec_tests.c > CMakeFiles/pastix_tests.dir/s_bvec_tests.c.i

test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/s_bvec_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_s -UPRECISION_p -UPRECISION_d -UPRECISION_c -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/s_bvec_tests.c -o CMakeFiles/pastix_tests.dir/s_bvec_tests.c.s

test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o


test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o: test/z_bvec_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_z -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_c -o CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_bvec_tests.c

test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/z_bvec_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_z -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_c -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_bvec_tests.c > CMakeFiles/pastix_tests.dir/z_bvec_tests.c.i

test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/z_bvec_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_z -UPRECISION_p -UPRECISION_s -UPRECISION_d -UPRECISION_c -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/z_bvec_tests.c -o CMakeFiles/pastix_tests.dir/z_bvec_tests.c.s

test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o


test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o: test/CMakeFiles/pastix_tests.dir/flags.make
test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o: test/d_bvec_tests.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -o CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o   -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_bvec_tests.c

test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pastix_tests.dir/d_bvec_tests.c.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_bvec_tests.c > CMakeFiles/pastix_tests.dir/d_bvec_tests.c.i

test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pastix_tests.dir/d_bvec_tests.c.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -DPRECISION_d -UPRECISION_p -UPRECISION_s -UPRECISION_c -UPRECISION_z -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/d_bvec_tests.c -o CMakeFiles/pastix_tests.dir/d_bvec_tests.c.s

test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.requires:

.PHONY : test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.requires

test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.provides: test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.requires
	$(MAKE) -f test/CMakeFiles/pastix_tests.dir/build.make test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.provides.build
.PHONY : test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.provides

test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.provides.build: test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o


# Object files for target pastix_tests
pastix_tests_OBJECTS = \
"CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o" \
"CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o" \
"CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o" \
"CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o" \
"CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o" \
"CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o" \
"CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o" \
"CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o"

# External object files for target pastix_tests
pastix_tests_EXTERNAL_OBJECTS =

test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/build.make
test/libpastix_tests.a: test/CMakeFiles/pastix_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Linking C static library libpastix_tests.a"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && $(CMAKE_COMMAND) -P CMakeFiles/pastix_tests.dir/cmake_clean_target.cmake
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pastix_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/pastix_tests.dir/build: test/libpastix_tests.a

.PHONY : test/CMakeFiles/pastix_tests.dir/build

test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/c_lowrank_tests.c.o.requires
test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/s_lowrank_tests.c.o.requires
test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/z_lowrank_tests.c.o.requires
test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/d_lowrank_tests.c.o.requires
test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/c_bvec_tests.c.o.requires
test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/s_bvec_tests.c.o.requires
test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/z_bvec_tests.c.o.requires
test/CMakeFiles/pastix_tests.dir/requires: test/CMakeFiles/pastix_tests.dir/d_bvec_tests.c.o.requires

.PHONY : test/CMakeFiles/pastix_tests.dir/requires

test/CMakeFiles/pastix_tests.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test && $(CMAKE_COMMAND) -P CMakeFiles/pastix_tests.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/pastix_tests.dir/clean

test/CMakeFiles/pastix_tests.dir/depend: test/c_lowrank_tests.c
test/CMakeFiles/pastix_tests.dir/depend: test/s_lowrank_tests.c
test/CMakeFiles/pastix_tests.dir/depend: test/d_lowrank_tests.c
test/CMakeFiles/pastix_tests.dir/depend: test/c_bvec_tests.c
test/CMakeFiles/pastix_tests.dir/depend: test/s_bvec_tests.c
test/CMakeFiles/pastix_tests.dir/depend: test/d_bvec_tests.c
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/test/CMakeFiles/pastix_tests.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/pastix_tests.dir/depend

