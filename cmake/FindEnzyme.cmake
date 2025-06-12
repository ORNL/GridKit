# 
#[[

Finds Enzyme Clang plugin

User may set:
- ENZYME_DIR

Author(s):
- Asher Mancinelli <ashermancinelli@gmail.com>
- Nicholson Koukpaizan <koukpaizannk@ornl.gov>

]]

find_package(Enzyme REQUIRED CONFIG 
             PATHS 
             ${ENZYME_DIR}
             ${ENZYME_DIR}/lib/cmake/Enzyme)
message(STATUS "Enzyme configuration found: ${Enzyme_CONFIG}")

find_library(ENZYME_LLVM_PLUGIN_LIBRARY
  NAMES
  LLVMEnzyme-${Enzyme_LLVM_VERSION_MAJOR}.so
  LLVMEnzyme-${Enzyme_LLVM_VERSION_MAJOR}.dylib
  LLVMEnzyme-${Enzyme_LLVM_VERSION_MAJOR}.dll
  PATHS
  ${ENZYME_DIR}
  ENV LD_LIBRARY_PATH 
  ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib
  REQUIRED)
message(STATUS "Enzyme LLVM plugin library: ${ENZYME_LLVM_PLUGIN_LIBRARY}")

find_library(ENZYME_CLANG_PLUGIN_LIBRARY
  NAMES
  ClangEnzyme-${Enzyme_LLVM_VERSION_MAJOR}.so
  ClangEnzyme-${Enzyme_LLVM_VERSION_MAJOR}.dylib
  ClangEnzyme-${Enzyme_LLVM_VERSION_MAJOR}.dll
  PATHS
  ${ENZYME_DIR}
  ENV LD_LIBRARY_PATH 
  ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib
  REQUIRED)
message(STATUS "Enzyme Clang plugin library: ${ENZYME_CLANG_PLUGIN_LIBRARY}")

find_program(GRIDKIT_LLVM_LINK llvm-link
             PATHS ${Enzyme_LLVM_BINARY_DIR}
             PATH_SUFFIXES
             bin
             REQUIRED)
message(STATUS "llvm-link: ${GRIDKIT_LLVM_LINK}")

find_program(GRIDKIT_OPT opt
             PATHS ${Enzyme_LLVM_BINARY_DIR}
             PATH_SUFFIXES
             bin
             REQUIRED)
message(STATUS "opt: ${GRIDKIT_OPT}")
