cmake_minimum_required(VERSION 3.16...3.27 FATAL_ERROR)

project(Sr90Plate)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

add_executable(sim sim.cc ${sources})
target_link_libraries(sim ${Geant4_LIBRARIES})

# Copy macros next to the executable so ./sim run.mac works from the build dir
file(GLOB macros ${PROJECT_SOURCE_DIR}/macros/*.mac)
file(COPY ${macros} DESTINATION ${PROJECT_BINARY_DIR})
