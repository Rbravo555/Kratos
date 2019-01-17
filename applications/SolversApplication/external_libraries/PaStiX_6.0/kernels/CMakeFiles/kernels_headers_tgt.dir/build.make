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

# Utility rule file for kernels_headers_tgt.

# Include the progress variables for this target.
include kernels/CMakeFiles/kernels_headers_tgt.dir/progress.make

kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_ccores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_scores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_zcores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_dcores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_ccuda.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_scuda.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_zcuda.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_dcuda.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_clrcores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_slrcores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_zlrcores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/pastix_dlrcores.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/c_nan_check.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/s_nan_check.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/z_nan_check.h
kernels/CMakeFiles/kernels_headers_tgt: kernels/d_nan_check.h


kernels/pastix_ccores.h: kernels/pastix_zcores.h
kernels/pastix_ccores.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_ccores.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating pastix_ccores.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_ccores.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zcores.h -p c -P ./ && chmod a-w pastix_ccores.h

kernels/pastix_scores.h: kernels/pastix_zcores.h
kernels/pastix_scores.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_scores.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating pastix_scores.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_scores.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zcores.h -p s -P ./ && chmod a-w pastix_scores.h

kernels/pastix_dcores.h: kernels/pastix_zcores.h
kernels/pastix_dcores.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_dcores.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating pastix_dcores.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_dcores.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zcores.h -p d -P ./ && chmod a-w pastix_dcores.h

kernels/pastix_ccuda.h: kernels/pastix_zcuda.h
kernels/pastix_ccuda.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_ccuda.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Generating pastix_ccuda.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_ccuda.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zcuda.h -p c -P ./ && chmod a-w pastix_ccuda.h

kernels/pastix_scuda.h: kernels/pastix_zcuda.h
kernels/pastix_scuda.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_scuda.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Generating pastix_scuda.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_scuda.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zcuda.h -p s -P ./ && chmod a-w pastix_scuda.h

kernels/pastix_dcuda.h: kernels/pastix_zcuda.h
kernels/pastix_dcuda.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_dcuda.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Generating pastix_dcuda.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_dcuda.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zcuda.h -p d -P ./ && chmod a-w pastix_dcuda.h

kernels/pastix_clrcores.h: kernels/pastix_zlrcores.h
kernels/pastix_clrcores.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_clrcores.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Generating pastix_clrcores.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_clrcores.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zlrcores.h -p c -P ./ && chmod a-w pastix_clrcores.h

kernels/pastix_slrcores.h: kernels/pastix_zlrcores.h
kernels/pastix_slrcores.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_slrcores.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Generating pastix_slrcores.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_slrcores.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zlrcores.h -p s -P ./ && chmod a-w pastix_slrcores.h

kernels/pastix_dlrcores.h: kernels/pastix_zlrcores.h
kernels/pastix_dlrcores.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/pastix_dlrcores.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Generating pastix_dlrcores.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f pastix_dlrcores.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/pastix_zlrcores.h -p d -P ./ && chmod a-w pastix_dlrcores.h

kernels/c_nan_check.h: kernels/z_nan_check.h
kernels/c_nan_check.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/c_nan_check.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Generating c_nan_check.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f c_nan_check.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/z_nan_check.h -p c -P ./ && chmod a-w c_nan_check.h

kernels/s_nan_check.h: kernels/z_nan_check.h
kernels/s_nan_check.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/s_nan_check.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Generating s_nan_check.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f s_nan_check.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/z_nan_check.h -p s -P ./ && chmod a-w s_nan_check.h

kernels/d_nan_check.h: kernels/z_nan_check.h
kernels/d_nan_check.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
kernels/d_nan_check.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Generating d_nan_check.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && /usr/bin/cmake -E remove -f d_nan_check.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/z_nan_check.h -p d -P ./ && chmod a-w d_nan_check.h

kernels_headers_tgt: kernels/CMakeFiles/kernels_headers_tgt
kernels_headers_tgt: kernels/pastix_ccores.h
kernels_headers_tgt: kernels/pastix_scores.h
kernels_headers_tgt: kernels/pastix_dcores.h
kernels_headers_tgt: kernels/pastix_ccuda.h
kernels_headers_tgt: kernels/pastix_scuda.h
kernels_headers_tgt: kernels/pastix_dcuda.h
kernels_headers_tgt: kernels/pastix_clrcores.h
kernels_headers_tgt: kernels/pastix_slrcores.h
kernels_headers_tgt: kernels/pastix_dlrcores.h
kernels_headers_tgt: kernels/c_nan_check.h
kernels_headers_tgt: kernels/s_nan_check.h
kernels_headers_tgt: kernels/d_nan_check.h
kernels_headers_tgt: kernels/CMakeFiles/kernels_headers_tgt.dir/build.make

.PHONY : kernels_headers_tgt

# Rule to build all files generated by this target.
kernels/CMakeFiles/kernels_headers_tgt.dir/build: kernels_headers_tgt

.PHONY : kernels/CMakeFiles/kernels_headers_tgt.dir/build

kernels/CMakeFiles/kernels_headers_tgt.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels && $(CMAKE_COMMAND) -P CMakeFiles/kernels_headers_tgt.dir/cmake_clean.cmake
.PHONY : kernels/CMakeFiles/kernels_headers_tgt.dir/clean

kernels/CMakeFiles/kernels_headers_tgt.dir/depend:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/kernels/CMakeFiles/kernels_headers_tgt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : kernels/CMakeFiles/kernels_headers_tgt.dir/depend

