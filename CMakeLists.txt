cmake_minimum_required(VERSION 3.16)
project(LexiPairBioRange CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# Slide 17: GCC -O3, single-threaded, high_resolution_clock.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
include_directories(include)

add_library(lexipair STATIC
  src/lexipair_index.cpp
  src/lsm_index.cpp
  src/fractional_cascade.cpp
)

foreach(app build_index static_query dynamic_update fc_benchmark scalability ablation vgp_compare validate)
  add_executable(${app} apps/${app}.cpp)
  target_link_libraries(${app} lexipair)
endforeach()
target_compile_definitions(scalability PRIVATE APP_SCALABILITY)
target_compile_definitions(ablation PRIVATE APP_ABLATION)
target_compile_definitions(vgp_compare PRIVATE APP_VGP)
