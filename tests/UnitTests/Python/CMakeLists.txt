find_package(Python3 COMPONENTS Interpreter)

if(Python3_FOUND)
  # Uses the venv built by the notebooks_pyenv target (see
  # application/notebooks/CMakeLists.txt) rather than the bare Python3
  # interpreter, since that venv is what actually has pytest and the
  # notebook packages' runtime dependencies installed. As with
  # `test_install` (see INSTALL.md), this target must be built before it
  # can run: `cmake --build . --target notebooks_pyenv` then `ctest`.
  if(WIN32)
    set(_nb_venv_python ${CMAKE_BINARY_DIR}/application/notebooks/notebooks-venv/Scripts/python.exe)
  else()
    set(_nb_venv_python ${CMAKE_BINARY_DIR}/application/notebooks/notebooks-venv/bin/python)
  endif()

  add_test(
    NAME PythonNotebookUtilsTest
    COMMAND ${_nb_venv_python} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(PythonNotebookUtilsTest PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/application/notebooks")
else()
  message(STATUS "Python3 interpreter not found; PythonNotebookUtilsTest unavailable")
endif()
