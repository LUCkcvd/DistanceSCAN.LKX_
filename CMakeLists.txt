cmake_minimum_required(VERSION 3.9)
project(Distance_SCAN_SIGMOD)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Add debug flags without -Werror
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Add OpenMP support
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
endif()

# Set output path
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

# Configure Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_NO_BOOST_CMAKE ON)

# Find Boost packages
find_package(Boost REQUIRED COMPONENTS 
    unit_test_framework 
    filesystem 
    system 
    date_time 
    program_options 
    serialization
)

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
    ${OpenMP_CXX_INCLUDE_DIRS}
    /usr/lib/gcc/aarch64-linux-gnu/14/include  # Add GCC include path for OpenMP
)

# Add compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Add source files
set(SOURCES
    main.cpp
    src/config.cpp
    src/graph.cpp
    src/mylib.cpp
    src/clusteralgo.cpp
    src/sketches.cpp
    src/fhq_treap.cpp
    src/binary_tree.cpp
)

# Add header files
set(HEADERS
    include/query.h
    include/clusteralgo.h
    include/sketches.h
    include/fhq_treap.h
    include/algos.h
    include/binary_tree.h
)

# Create executable
add_executable(Distance_SCAN_SIGMOD ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(Distance_SCAN_SIGMOD
    PRIVATE
    ${Boost_LIBRARIES}
    OpenMP::OpenMP_CXX
    gomp
    pthread
)

# Enable debugging symbols
set_target_properties(Distance_SCAN_SIGMOD
    PROPERTIES
    COMPILE_FLAGS "-g"
    LINK_FLAGS "-g"
)
