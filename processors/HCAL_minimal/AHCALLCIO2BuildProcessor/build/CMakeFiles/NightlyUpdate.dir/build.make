# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/CMake/3.15.5/bin/cmake

# The command to remove a file.
RM = /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/CMake/3.15.5/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor/build

# Utility rule file for NightlyUpdate.

# Include the progress variables for this target.
include CMakeFiles/NightlyUpdate.dir/progress.make

CMakeFiles/NightlyUpdate:
	/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/CMake/3.15.5/bin/ctest -D NightlyUpdate

NightlyUpdate: CMakeFiles/NightlyUpdate
NightlyUpdate: CMakeFiles/NightlyUpdate.dir/build.make

.PHONY : NightlyUpdate

# Rule to build all files generated by this target.
CMakeFiles/NightlyUpdate.dir/build: NightlyUpdate

.PHONY : CMakeFiles/NightlyUpdate.dir/build

CMakeFiles/NightlyUpdate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NightlyUpdate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NightlyUpdate.dir/clean

CMakeFiles/NightlyUpdate.dir/depend:
	cd /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor/build /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor/build /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/HCAL_minimal/AHCALLCIO2BuildProcessor/build/CMakeFiles/NightlyUpdate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NightlyUpdate.dir/depend

