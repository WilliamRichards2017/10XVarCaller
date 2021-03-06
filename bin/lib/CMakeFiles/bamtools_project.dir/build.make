# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /uufs/chpc.utah.edu/common/home/u0401321/TenX

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin

# Utility rule file for bamtools_project.

# Include the progress variables for this target.
include lib/CMakeFiles/bamtools_project.dir/progress.make

lib/CMakeFiles/bamtools_project: lib/CMakeFiles/bamtools_project-complete

lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-install
lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-mkdir
lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-download
lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-update
lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-patch
lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure
lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-build
lib/CMakeFiles/bamtools_project-complete: externals/bamtools/src/bamtools_project-stamp/bamtools_project-install
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Completed 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E make_directory /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib/CMakeFiles
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib/CMakeFiles/bamtools_project-complete
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-done

externals/bamtools/src/bamtools_project-stamp/bamtools_project-install: externals/bamtools/src/bamtools_project-stamp/bamtools_project-build
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No install step for 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-build && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-install

externals/bamtools/src/bamtools_project-stamp/bamtools_project-mkdir:
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Creating directories for 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E make_directory /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E make_directory /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-build
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E make_directory /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E make_directory /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/tmp
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E make_directory /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E make_directory /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-mkdir

externals/bamtools/src/bamtools_project-stamp/bamtools_project-download: externals/bamtools/src/bamtools_project-stamp/bamtools_project-gitinfo.txt
externals/bamtools/src/bamtools_project-stamp/bamtools_project-download: externals/bamtools/src/bamtools_project-stamp/bamtools_project-mkdir
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing download step (git clone) for 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src && /usr/bin/cmake -P /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/tmp/bamtools_project-gitclone.cmake
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-download

externals/bamtools/src/bamtools_project-stamp/bamtools_project-update: externals/bamtools/src/bamtools_project-stamp/bamtools_project-download
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No update step for 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-update

externals/bamtools/src/bamtools_project-stamp/bamtools_project-patch: externals/bamtools/src/bamtools_project-stamp/bamtools_project-download
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No patch step for 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-patch

externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure: externals/bamtools/tmp/bamtools_project-cfgcmd.txt
externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure: externals/bamtools/tmp/bamtools_project-cache.cmake
externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure: externals/bamtools/src/bamtools_project-stamp/bamtools_project-update
externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure: externals/bamtools/src/bamtools_project-stamp/bamtools_project-patch
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing configure step for 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-build && /usr/bin/cmake -C/uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/tmp/bamtools_project-cache.cmake "-GUnix Makefiles" /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-build && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure

externals/bamtools/src/bamtools_project-stamp/bamtools_project-build: externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure
	$(CMAKE_COMMAND) -E cmake_progress_report /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing build step for 'bamtools_project'"
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-build && $(MAKE)
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-build && /usr/bin/cmake -E touch /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/externals/bamtools/src/bamtools_project-stamp/bamtools_project-build

bamtools_project: lib/CMakeFiles/bamtools_project
bamtools_project: lib/CMakeFiles/bamtools_project-complete
bamtools_project: externals/bamtools/src/bamtools_project-stamp/bamtools_project-install
bamtools_project: externals/bamtools/src/bamtools_project-stamp/bamtools_project-mkdir
bamtools_project: externals/bamtools/src/bamtools_project-stamp/bamtools_project-download
bamtools_project: externals/bamtools/src/bamtools_project-stamp/bamtools_project-update
bamtools_project: externals/bamtools/src/bamtools_project-stamp/bamtools_project-patch
bamtools_project: externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure
bamtools_project: externals/bamtools/src/bamtools_project-stamp/bamtools_project-build
bamtools_project: lib/CMakeFiles/bamtools_project.dir/build.make
.PHONY : bamtools_project

# Rule to build all files generated by this target.
lib/CMakeFiles/bamtools_project.dir/build: bamtools_project
.PHONY : lib/CMakeFiles/bamtools_project.dir/build

lib/CMakeFiles/bamtools_project.dir/clean:
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib && $(CMAKE_COMMAND) -P CMakeFiles/bamtools_project.dir/cmake_clean.cmake
.PHONY : lib/CMakeFiles/bamtools_project.dir/clean

lib/CMakeFiles/bamtools_project.dir/depend:
	cd /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /uufs/chpc.utah.edu/common/home/u0401321/TenX /uufs/chpc.utah.edu/common/home/u0401321/TenX/lib /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib /uufs/chpc.utah.edu/common/home/u0401321/TenX/bin/lib/CMakeFiles/bamtools_project.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/CMakeFiles/bamtools_project.dir/depend

