# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = C:\Users\bruno\AppData\Local\JetBrains\Toolbox\apps\CLion\ch-0\202.6948.80\bin\cmake\win\bin\cmake.exe

# The command to remove a file.
RM = C:\Users\bruno\AppData\Local\JetBrains\Toolbox\apps\CLion\ch-0\202.6948.80\bin\cmake\win\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\bruno\Documents\Github\CorridorAllocation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\bruno\Documents\Github\CorridorAllocation\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/OrgSalas.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/OrgSalas.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/OrgSalas.dir/flags.make

CMakeFiles/OrgSalas.dir/main.cpp.obj: CMakeFiles/OrgSalas.dir/flags.make
CMakeFiles/OrgSalas.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\bruno\Documents\Github\CorridorAllocation\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/OrgSalas.dir/main.cpp.obj"
	C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\OrgSalas.dir\main.cpp.obj -c C:\Users\bruno\Documents\Github\CorridorAllocation\main.cpp

CMakeFiles/OrgSalas.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OrgSalas.dir/main.cpp.i"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\bruno\Documents\Github\CorridorAllocation\main.cpp > CMakeFiles\OrgSalas.dir\main.cpp.i

CMakeFiles/OrgSalas.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OrgSalas.dir/main.cpp.s"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\bruno\Documents\Github\CorridorAllocation\main.cpp -o CMakeFiles\OrgSalas.dir\main.cpp.s

# Object files for target OrgSalas
OrgSalas_OBJECTS = \
"CMakeFiles/OrgSalas.dir/main.cpp.obj"

# External object files for target OrgSalas
OrgSalas_EXTERNAL_OBJECTS =

OrgSalas.exe: CMakeFiles/OrgSalas.dir/main.cpp.obj
OrgSalas.exe: CMakeFiles/OrgSalas.dir/build.make
OrgSalas.exe: CMakeFiles/OrgSalas.dir/linklibs.rsp
OrgSalas.exe: CMakeFiles/OrgSalas.dir/objects1.rsp
OrgSalas.exe: CMakeFiles/OrgSalas.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\bruno\Documents\Github\CorridorAllocation\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable OrgSalas.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\OrgSalas.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/OrgSalas.dir/build: OrgSalas.exe

.PHONY : CMakeFiles/OrgSalas.dir/build

CMakeFiles/OrgSalas.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\OrgSalas.dir\cmake_clean.cmake
.PHONY : CMakeFiles/OrgSalas.dir/clean

CMakeFiles/OrgSalas.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\bruno\Documents\Github\CorridorAllocation C:\Users\bruno\Documents\Github\CorridorAllocation C:\Users\bruno\Documents\Github\CorridorAllocation\cmake-build-debug C:\Users\bruno\Documents\Github\CorridorAllocation\cmake-build-debug C:\Users\bruno\Documents\Github\CorridorAllocation\cmake-build-debug\CMakeFiles\OrgSalas.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/OrgSalas.dir/depend

