# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.25.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.25.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build

# Include any dependencies generated for this target.
include CMakeFiles/rdecay01.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/rdecay01.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/rdecay01.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/rdecay01.dir/flags.make

CMakeFiles/rdecay01.dir/rdecay01.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/rdecay01.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/rdecay01.cc
CMakeFiles/rdecay01.dir/rdecay01.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/rdecay01.dir/rdecay01.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/rdecay01.cc.o -MF CMakeFiles/rdecay01.dir/rdecay01.cc.o.d -o CMakeFiles/rdecay01.dir/rdecay01.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/rdecay01.cc

CMakeFiles/rdecay01.dir/rdecay01.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/rdecay01.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/rdecay01.cc > CMakeFiles/rdecay01.dir/rdecay01.cc.i

CMakeFiles/rdecay01.dir/rdecay01.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/rdecay01.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/rdecay01.cc -o CMakeFiles/rdecay01.dir/rdecay01.cc.s

CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/ActionInitialization.cc
CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o -MF CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o.d -o CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/ActionInitialization.cc

CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/ActionInitialization.cc > CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.i

CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/ActionInitialization.cc -o CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.s

CMakeFiles/rdecay01.dir/src/Detector.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/Detector.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Detector.cc
CMakeFiles/rdecay01.dir/src/Detector.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/rdecay01.dir/src/Detector.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/Detector.cc.o -MF CMakeFiles/rdecay01.dir/src/Detector.cc.o.d -o CMakeFiles/rdecay01.dir/src/Detector.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Detector.cc

CMakeFiles/rdecay01.dir/src/Detector.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/Detector.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Detector.cc > CMakeFiles/rdecay01.dir/src/Detector.cc.i

CMakeFiles/rdecay01.dir/src/Detector.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/Detector.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Detector.cc -o CMakeFiles/rdecay01.dir/src/Detector.cc.s

CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/DetectorConstruction.cc
CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o -MF CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o.d -o CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/DetectorConstruction.cc

CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/DetectorConstruction.cc > CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.i

CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/DetectorConstruction.cc -o CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.s

CMakeFiles/rdecay01.dir/src/EventAction.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/EventAction.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/EventAction.cc
CMakeFiles/rdecay01.dir/src/EventAction.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/rdecay01.dir/src/EventAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/EventAction.cc.o -MF CMakeFiles/rdecay01.dir/src/EventAction.cc.o.d -o CMakeFiles/rdecay01.dir/src/EventAction.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/EventAction.cc

CMakeFiles/rdecay01.dir/src/EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/EventAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/EventAction.cc > CMakeFiles/rdecay01.dir/src/EventAction.cc.i

CMakeFiles/rdecay01.dir/src/EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/EventAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/EventAction.cc -o CMakeFiles/rdecay01.dir/src/EventAction.cc.s

CMakeFiles/rdecay01.dir/src/HistoManager.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/HistoManager.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/HistoManager.cc
CMakeFiles/rdecay01.dir/src/HistoManager.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/rdecay01.dir/src/HistoManager.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/HistoManager.cc.o -MF CMakeFiles/rdecay01.dir/src/HistoManager.cc.o.d -o CMakeFiles/rdecay01.dir/src/HistoManager.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/HistoManager.cc

CMakeFiles/rdecay01.dir/src/HistoManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/HistoManager.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/HistoManager.cc > CMakeFiles/rdecay01.dir/src/HistoManager.cc.i

CMakeFiles/rdecay01.dir/src/HistoManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/HistoManager.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/HistoManager.cc -o CMakeFiles/rdecay01.dir/src/HistoManager.cc.s

CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PhysicsList.cc
CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o -MF CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o.d -o CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PhysicsList.cc

CMakeFiles/rdecay01.dir/src/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/PhysicsList.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PhysicsList.cc > CMakeFiles/rdecay01.dir/src/PhysicsList.cc.i

CMakeFiles/rdecay01.dir/src/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/PhysicsList.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PhysicsList.cc -o CMakeFiles/rdecay01.dir/src/PhysicsList.cc.s

CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PrimaryGeneratorAction.cc
CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o -MF CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o.d -o CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PrimaryGeneratorAction.cc

CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PrimaryGeneratorAction.cc > CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.i

CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/PrimaryGeneratorAction.cc -o CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.s

CMakeFiles/rdecay01.dir/src/Run.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/Run.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Run.cc
CMakeFiles/rdecay01.dir/src/Run.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/rdecay01.dir/src/Run.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/Run.cc.o -MF CMakeFiles/rdecay01.dir/src/Run.cc.o.d -o CMakeFiles/rdecay01.dir/src/Run.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Run.cc

CMakeFiles/rdecay01.dir/src/Run.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/Run.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Run.cc > CMakeFiles/rdecay01.dir/src/Run.cc.i

CMakeFiles/rdecay01.dir/src/Run.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/Run.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/Run.cc -o CMakeFiles/rdecay01.dir/src/Run.cc.s

CMakeFiles/rdecay01.dir/src/RunAction.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/RunAction.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/RunAction.cc
CMakeFiles/rdecay01.dir/src/RunAction.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/rdecay01.dir/src/RunAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/RunAction.cc.o -MF CMakeFiles/rdecay01.dir/src/RunAction.cc.o.d -o CMakeFiles/rdecay01.dir/src/RunAction.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/RunAction.cc

CMakeFiles/rdecay01.dir/src/RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/RunAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/RunAction.cc > CMakeFiles/rdecay01.dir/src/RunAction.cc.i

CMakeFiles/rdecay01.dir/src/RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/RunAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/RunAction.cc -o CMakeFiles/rdecay01.dir/src/RunAction.cc.s

CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingAction.cc
CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o -MF CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o.d -o CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingAction.cc

CMakeFiles/rdecay01.dir/src/TrackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/TrackingAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingAction.cc > CMakeFiles/rdecay01.dir/src/TrackingAction.cc.i

CMakeFiles/rdecay01.dir/src/TrackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/TrackingAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingAction.cc -o CMakeFiles/rdecay01.dir/src/TrackingAction.cc.s

CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o: CMakeFiles/rdecay01.dir/flags.make
CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o: /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingMessenger.cc
CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o: CMakeFiles/rdecay01.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o -MF CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o.d -o CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o -c /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingMessenger.cc

CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingMessenger.cc > CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.i

CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/src/TrackingMessenger.cc -o CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.s

# Object files for target rdecay01
rdecay01_OBJECTS = \
"CMakeFiles/rdecay01.dir/rdecay01.cc.o" \
"CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/rdecay01.dir/src/Detector.cc.o" \
"CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/rdecay01.dir/src/EventAction.cc.o" \
"CMakeFiles/rdecay01.dir/src/HistoManager.cc.o" \
"CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o" \
"CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/rdecay01.dir/src/Run.cc.o" \
"CMakeFiles/rdecay01.dir/src/RunAction.cc.o" \
"CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o" \
"CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o"

# External object files for target rdecay01
rdecay01_EXTERNAL_OBJECTS =

rdecay01: CMakeFiles/rdecay01.dir/rdecay01.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/ActionInitialization.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/Detector.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/DetectorConstruction.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/EventAction.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/HistoManager.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/PhysicsList.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/PrimaryGeneratorAction.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/Run.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/RunAction.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/TrackingAction.cc.o
rdecay01: CMakeFiles/rdecay01.dir/src/TrackingMessenger.cc.o
rdecay01: CMakeFiles/rdecay01.dir/build.make
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4Tree.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4FR.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4GMocren.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4visHepRep.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4RayTracer.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4VRML.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4OpenGL.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4gl2ps.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4visQt3D.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4ToolsSG.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4vis_management.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4modeling.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4interfaces.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4persistency.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4error_propagation.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4readout.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4physicslists.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4run.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4event.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4tracking.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4parmodels.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4processes.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4digits_hits.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4track.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4particles.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4geometry.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4materials.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4graphics_reps.dylib
rdecay01: /usr/local/opt/qt5/lib/Qt3DExtras.framework/Qt3DExtras
rdecay01: /usr/local/opt/qt5/lib/Qt3DInput.framework/Qt3DInput
rdecay01: /usr/local/opt/qt5/lib/QtGamepad.framework/QtGamepad
rdecay01: /usr/local/opt/qt5/lib/Qt3DLogic.framework/Qt3DLogic
rdecay01: /usr/local/opt/qt5/lib/Qt3DRender.framework/Qt3DRender
rdecay01: /usr/local/opt/qt5/lib/Qt3DCore.framework/Qt3DCore
rdecay01: /usr/local/opt/qt5/lib/QtNetwork.framework/QtNetwork
rdecay01: /Library/Developer/CommandLineTools/SDKs/MacOSX11.0.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd
rdecay01: /usr/local/opt/qt5/lib/QtOpenGL.framework/QtOpenGL
rdecay01: /usr/local/opt/qt5/lib/QtPrintSupport.framework/QtPrintSupport
rdecay01: /usr/local/opt/qt5/lib/QtWidgets.framework/QtWidgets
rdecay01: /usr/local/opt/qt5/lib/QtGui.framework/QtGui
rdecay01: /usr/local/opt/qt5/lib/QtCore.framework/QtCore
rdecay01: /usr/local/lib/libxerces-c.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4analysis.dylib
rdecay01: /Library/Developer/CommandLineTools/SDKs/MacOSX11.0.sdk/usr/lib/libexpat.tbd
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4zlib.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4intercoms.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4global.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4clhep.dylib
rdecay01: /Users/samueletorelli/Documents/geant4/build/BuildProducts/lib/libG4ptl.2.3.3.dylib
rdecay01: CMakeFiles/rdecay01.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable rdecay01"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rdecay01.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/rdecay01.dir/build: rdecay01
.PHONY : CMakeFiles/rdecay01.dir/build

CMakeFiles/rdecay01.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/rdecay01.dir/cmake_clean.cmake
.PHONY : CMakeFiles/rdecay01.dir/clean

CMakeFiles/rdecay01.dir/depend:
	cd /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build /Users/samueletorelli/Documents/Dottorato/CYGNO_30/CYGNO_30_Sim/build/CMakeFiles/rdecay01.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/rdecay01.dir/depend
