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

# Utility rule file for spm_headers_tgt.

# Include the progress variables for this target.
include spm/CMakeFiles/spm_headers_tgt.dir/progress.make

spm/CMakeFiles/spm_headers_tgt: spm/include/p_spm.h
spm/CMakeFiles/spm_headers_tgt: spm/include/s_spm.h
spm/CMakeFiles/spm_headers_tgt: spm/include/c_spm.h
spm/CMakeFiles/spm_headers_tgt: spm/include/d_spm.h
spm/CMakeFiles/spm_headers_tgt: spm/include/z_spm.h
spm/CMakeFiles/spm_headers_tgt: spm/include/spm.h
spm/CMakeFiles/spm_headers_tgt: spm/src/spm_drivers.h


spm/include/p_spm.h: spm/src/z_spm.h
spm/include/p_spm.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
spm/include/p_spm.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating include/p_spm.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm && /usr/bin/cmake -E remove -f include/p_spm.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/src/z_spm.h -p p -P include && chmod a-w include/p_spm.h

spm/include/s_spm.h: spm/src/z_spm.h
spm/include/s_spm.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
spm/include/s_spm.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating include/s_spm.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm && /usr/bin/cmake -E remove -f include/s_spm.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/src/z_spm.h -p s -P include && chmod a-w include/s_spm.h

spm/include/c_spm.h: spm/src/z_spm.h
spm/include/c_spm.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
spm/include/c_spm.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating include/c_spm.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm && /usr/bin/cmake -E remove -f include/c_spm.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/src/z_spm.h -p c -P include && chmod a-w include/c_spm.h

spm/include/d_spm.h: spm/src/z_spm.h
spm/include/d_spm.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
spm/include/d_spm.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Generating include/d_spm.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm && /usr/bin/cmake -E remove -f include/d_spm.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/src/z_spm.h -p d -P include && chmod a-w include/d_spm.h

spm/include/z_spm.h: spm/src/z_spm.h
spm/include/z_spm.h: cmake_modules/morse_cmake/modules/precision_generator/codegen.py
spm/include/z_spm.h: cmake_modules/morse_cmake/modules/precision_generator/subs.py
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Generating include/z_spm.h"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm && /usr/bin/cmake -E remove -f include/z_spm.h && /usr/bin/python2 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/cmake_modules/morse_cmake/modules/precision_generator/codegen.py -f /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/src/z_spm.h -p z -P include && chmod a-w include/z_spm.h

spm_headers_tgt: spm/CMakeFiles/spm_headers_tgt
spm_headers_tgt: spm/include/p_spm.h
spm_headers_tgt: spm/include/s_spm.h
spm_headers_tgt: spm/include/c_spm.h
spm_headers_tgt: spm/include/d_spm.h
spm_headers_tgt: spm/include/z_spm.h
spm_headers_tgt: spm/CMakeFiles/spm_headers_tgt.dir/build.make

.PHONY : spm_headers_tgt

# Rule to build all files generated by this target.
spm/CMakeFiles/spm_headers_tgt.dir/build: spm_headers_tgt

.PHONY : spm/CMakeFiles/spm_headers_tgt.dir/build

spm/CMakeFiles/spm_headers_tgt.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm && $(CMAKE_COMMAND) -P CMakeFiles/spm_headers_tgt.dir/cmake_clean.cmake
.PHONY : spm/CMakeFiles/spm_headers_tgt.dir/clean

spm/CMakeFiles/spm_headers_tgt.dir/depend:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/CMakeFiles/spm_headers_tgt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : spm/CMakeFiles/spm_headers_tgt.dir/depend

