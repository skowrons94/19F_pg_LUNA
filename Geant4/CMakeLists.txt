#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4_BGO)

set(CMAKE_CXX_STANDARD 14)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB_RECURSE sources RELATIVE ${PROJECT_SOURCE_DIR} *.cc)

file(GLOB_RECURSE headers RELATIVE ${PROJECT_SOURCE_DIR} *.hh)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(G4_BGO ${sources} ${headers})
target_link_libraries(G4_BGO ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory
#
file(COPY ${PROJECT_SOURCE_DIR}/mac DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS G4_BGO DESTINATION bin)
