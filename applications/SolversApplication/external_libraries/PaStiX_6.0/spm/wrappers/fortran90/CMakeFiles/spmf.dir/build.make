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
include spm/wrappers/fortran90/CMakeFiles/spmf.dir/depend.make

# Include the progress variables for this target.
include spm/wrappers/fortran90/CMakeFiles/spmf.dir/progress.make

# Include the compile flags for this target's objects.
include spm/wrappers/fortran90/CMakeFiles/spmf.dir/flags.make

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o: spm/wrappers/fortran90/CMakeFiles/spmf.dir/flags.make
spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o: spm/wrappers/fortran90/src/spm_enums.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && /usr/bin/gfortran $(Fortran_DEFINES) -DSPM_INT_KIND=c_int32_t $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/src/spm_enums.F90 -o CMakeFiles/spmf.dir/src/spm_enums.F90.o

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/spmf.dir/src/spm_enums.F90.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && /usr/bin/gfortran $(Fortran_DEFINES) -DSPM_INT_KIND=c_int32_t $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/src/spm_enums.F90 > CMakeFiles/spmf.dir/src/spm_enums.F90.i

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/spmf.dir/src/spm_enums.F90.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && /usr/bin/gfortran $(Fortran_DEFINES) -DSPM_INT_KIND=c_int32_t $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/src/spm_enums.F90 -o CMakeFiles/spmf.dir/src/spm_enums.F90.s

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.requires:

.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.requires

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.provides: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.requires
	$(MAKE) -f spm/wrappers/fortran90/CMakeFiles/spmf.dir/build.make spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.provides.build
.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.provides

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.provides.build: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o


spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o: spm/wrappers/fortran90/CMakeFiles/spmf.dir/flags.make
spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o: spm/wrappers/fortran90/src/spmf.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/src/spmf.f90 -o CMakeFiles/spmf.dir/src/spmf.f90.o

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/spmf.dir/src/spmf.f90.i"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/src/spmf.f90 > CMakeFiles/spmf.dir/src/spmf.f90.i

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/spmf.dir/src/spmf.f90.s"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && /usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/src/spmf.f90 -o CMakeFiles/spmf.dir/src/spmf.f90.s

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.requires:

.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.requires

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.provides: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.requires
	$(MAKE) -f spm/wrappers/fortran90/CMakeFiles/spmf.dir/build.make spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.provides.build
.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.provides

spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.provides.build: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o


# Object files for target spmf
spmf_OBJECTS = \
"CMakeFiles/spmf.dir/src/spm_enums.F90.o" \
"CMakeFiles/spmf.dir/src/spmf.f90.o"

# External object files for target spmf
spmf_EXTERNAL_OBJECTS =

spm/wrappers/fortran90/libspmf.a: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o
spm/wrappers/fortran90/libspmf.a: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o
spm/wrappers/fortran90/libspmf.a: spm/wrappers/fortran90/CMakeFiles/spmf.dir/build.make
spm/wrappers/fortran90/libspmf.a: spm/wrappers/fortran90/CMakeFiles/spmf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran static library libspmf.a"
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && $(CMAKE_COMMAND) -P CMakeFiles/spmf.dir/cmake_clean_target.cmake
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/spmf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
spm/wrappers/fortran90/CMakeFiles/spmf.dir/build: spm/wrappers/fortran90/libspmf.a

.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/build

spm/wrappers/fortran90/CMakeFiles/spmf.dir/requires: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spm_enums.F90.o.requires
spm/wrappers/fortran90/CMakeFiles/spmf.dir/requires: spm/wrappers/fortran90/CMakeFiles/spmf.dir/src/spmf.f90.o.requires

.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/requires

spm/wrappers/fortran90/CMakeFiles/spmf.dir/clean:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 && $(CMAKE_COMMAND) -P CMakeFiles/spmf.dir/cmake_clean.cmake
.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/clean

spm/wrappers/fortran90/CMakeFiles/spmf.dir/depend:
	cd /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90 /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/spm/wrappers/fortran90/CMakeFiles/spmf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : spm/wrappers/fortran90/CMakeFiles/spmf.dir/depend

