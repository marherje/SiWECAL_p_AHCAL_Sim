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
CMAKE_SOURCE_DIR = /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build

# Utility rule file for NightlyCoverage.

# Include the progress variables for this target.
include CMakeFiles/NightlyCoverage.dir/progress.make

CMakeFiles/NightlyCoverage:
	/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/CMake/3.15.5/bin/ctest -D NightlyCoverage

NightlyCoverage: CMakeFiles/NightlyCoverage
NightlyCoverage: CMakeFiles/NightlyCoverage.dir/build.make

.PHONY : NightlyCoverage

# Rule to build all files generated by this target.
CMakeFiles/NightlyCoverage.dir/build: NightlyCoverage

.PHONY : CMakeFiles/NightlyCoverage.dir/build

CMakeFiles/NightlyCoverage.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NightlyCoverage.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NightlyCoverage.dir/clean

CMakeFiles/NightlyCoverage.dir/depend:
	cd /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build /lhome/ific/m/marherje/SiWECAL_p_AHCAL_Sim/processors/ECAL/LCIO2BuildProcessor/build/CMakeFiles/NightlyCoverage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NightlyCoverage.dir/depend

