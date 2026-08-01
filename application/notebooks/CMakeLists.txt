option(GridKit_ENABLE_NOTEBOOKS_PYENV
       "Create a venv and install application/notebooks Python dependencies" OFF)

find_package(Python3 COMPONENTS Interpreter)

if(Python3_FOUND)
  set(_nb_venv ${CMAKE_CURRENT_BINARY_DIR}/notebooks-venv)

  add_custom_target(
    notebooks_pyenv
    COMMAND ${Python3_EXECUTABLE} -m venv ${_nb_venv}
    COMMAND ${_nb_venv}/bin/pip install --upgrade pip
    COMMAND ${_nb_venv}/bin/pip install -e ${CMAKE_CURRENT_SOURCE_DIR}[notebook,test]
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Creating venv and installing application/notebooks Python dependencies"
    VERBATIM)

  if(GridKit_ENABLE_NOTEBOOKS_PYENV)
    add_custom_target(notebooks_pyenv_all ALL DEPENDS notebooks_pyenv)
  endif()
else()
  message(STATUS "Python3 interpreter not found; notebooks_pyenv target unavailable")
endif()
