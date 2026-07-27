# Install script for directory: /home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/nopt/nrel/apps/cpu_stack/software/gridkit/llvm/16.0.6/bin/llvm-objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic"
         RPATH "/usr/local/lib:/nopt/nrel/apps/cpu_stack/software/gridkit/sundials/develop/lib64:/nopt/nrel/apps/cpu_stack/software/gridkit/suitesparse/7.8.2/lib64")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic" TYPE EXECUTABLE FILES "/home/isatkaus/gridkit/build/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic"
         OLD_RPATH "/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics:/home/isatkaus/gridkit/build/GridKit/Solver/Dynamic:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Branch:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Bus:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/BusFault:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/BusToSignalAdapter:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Exciter/IEEET1:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Governor/Tgov1:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Load/LoadZ:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Load/LoadZIP:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SignalSource:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Stabilizer/IEEEST:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SignalNode:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical:/nopt/nrel/apps/cpu_stack/software/gridkit/sundials/develop/lib64:/nopt/nrel/apps/cpu_stack/software/gridkit/suitesparse/7.8.2/lib64:/home/isatkaus/gridkit/build/GridKit/LinearAlgebra/SparseMatrix:/home/isatkaus/gridkit/build/GridKit/Utilities/Logger:/home/isatkaus/gridkit/build/GridKit/MemoryUtilities/cpu:"
         NEW_RPATH "/usr/local/lib:/nopt/nrel/apps/cpu_stack/software/gridkit/sundials/develop/lib64:/nopt/nrel/apps/cpu_stack/software/gridkit/suitesparse/7.8.2/lib64")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/nopt/nrel/apps/cpu_stack/software/gridkit/llvm/16.0.6/bin/llvm-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson"
         RPATH "/usr/local/lib:/nopt/nrel/apps/cpu_stack/software/gridkit/sundials/develop/lib64:/nopt/nrel/apps/cpu_stack/software/gridkit/suitesparse/7.8.2/lib64")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic" TYPE EXECUTABLE FILES "/home/isatkaus/gridkit/build/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson"
         OLD_RPATH "/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics:/home/isatkaus/gridkit/build/GridKit/Solver/Dynamic:/home/isatkaus/gridkit/build/GridKit/Testing:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Branch:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Bus:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/BusFault:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/BusToSignalAdapter:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Exciter/IEEET1:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Governor/Tgov1:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Load/LoadZ:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Load/LoadZIP:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SignalSource:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/Stabilizer/IEEEST:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SignalNode:/home/isatkaus/gridkit/build/GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical:/nopt/nrel/apps/cpu_stack/software/gridkit/sundials/develop/lib64:/nopt/nrel/apps/cpu_stack/software/gridkit/suitesparse/7.8.2/lib64:/home/isatkaus/gridkit/build/GridKit/LinearAlgebra/SparseMatrix:/home/isatkaus/gridkit/build/GridKit/MemoryUtilities/cpu:/home/isatkaus/gridkit/build/GridKit/Utilities/CliArgs:/home/isatkaus/gridkit/build/GridKit/Utilities/Logger:"
         NEW_RPATH "/usr/local/lib:/nopt/nrel/apps/cpu_stack/software/gridkit/sundials/develop/lib64:/nopt/nrel/apps/cpu_stack/software/gridkit/suitesparse/7.8.2/lib64")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/nopt/nrel/apps/cpu_stack/software/gridkit/llvm/16.0.6/bin/llvm-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasicJson")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic" TYPE FILE FILES "/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic.case.json")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic" TYPE FILE FILES "/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic.ref.csv")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic" TYPE FILE FILES "/home/isatkaus/gridkit/examples/PhasorDynamics/Tiny/ThreeBus/Basic/ThreeBusBasic.solver.json")
endif()

