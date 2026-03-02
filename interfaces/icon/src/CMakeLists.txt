cmake_minimum_required(VERSION 3.20)
project(p3plugin LANGUAGES C Fortran)

set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake")
message(STATUS "Looking for FindNetCDF.cmake in: ${CMAKE_MODULE_PATH}")

# use this for horeka module lib/netcdf-fortran/4.6: 
set(NF_CONFIG_BIN "/software/all/lib/netcdf-fortran/4.6_intel_2025.1_llvm_openmpi_5.0/bin")

# or instead set to whatever lies in $NETCDF_FORTRAN_BIN_DIR or $NETCDF_FORTRAN_ROOT/bin (example):
#set(NF_CONFIG_BIN "/software/all/lib/netcdf-fortran/4.6_intel_2025.1_llvm_impi_2021.11/bin")
#set(NF_CONFIG_BIN "${$NETCDF_FORTRAN_BIN_DIR}")  <-- sth like this does not work, cmake needs hard-coded paths

# compile the plugin with OpenMP for hybrid OpenMP-MPI parallelization:
#find_package(OpenMP REQUIRED COMPONENTS Fortran)

find_package(NetCDF REQUIRED COMPONENTS Fortran)
find_package(ComIn REQUIRED)

add_library(p3plugin_types SHARED p3plugin_types.F90)
add_library(p3plugin_global_vars SHARED p3plugin_global_vars.F90)
add_library(p3plugin_utils SHARED p3plugin_utils.F90)

add_library(p3plugin_tracer_init SHARED p3plugin_tracer_init.F90)
add_library(p3plugin_tracer_lbc SHARED p3plugin_tracer_lbc.F90)
add_library(p3plugin_add_tendencies SHARED p3plugin_add_tendencies.F90)
add_library(p3plugin_main_routines SHARED p3plugin_main_routines.F90)

add_library(p3_main SHARED p3_main.F90)
add_library(${PROJECT_NAME} MODULE p3plugin.F90)


# compiler options without runtime optimizations and debugging options:
target_compile_options(${PROJECT_NAME} PRIVATE
  -g -traceback -debug -O0 -fma -fp-model precise -ftz -qno-opt-dynamic-align)
target_compile_options(p3_main PRIVATE
  -g -traceback -debug -O0 -fma -fp-model precise -ftz -qno-opt-dynamic-align)

# compiler options with fast runtime optimizations without debugging:
#target_compile_options(${PROJECT_NAME} PRIVATE
#  -O2 -xICELAKE-SERVER -fma -fp-model fast -ftz -qno-opt-dynamic-align)
#target_compile_options(p3_main PRIVATE
#  -O2 -xICELAKE-SERVER -fma -fp-model fast -ftz -qno-opt-dynamic-align)


# compile the plugin with OpenMP for hybrid OpenMP-MPI parallelization:
#target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)

target_link_libraries(p3plugin_types PRIVATE ComIn::ComIn)
target_link_libraries(p3plugin_global_vars PRIVATE p3plugin_types ComIn::ComIn)
target_link_libraries(p3plugin_utils PRIVATE ComIn::ComIn)

target_link_libraries(p3plugin_tracer_init PRIVATE p3_main p3plugin_global_vars p3plugin_types ComIn::ComIn NetCDF::NetCDF_Fortran)
target_link_libraries(p3plugin_tracer_lbc PRIVATE p3plugin_global_vars p3plugin_types ComIn::ComIn)
target_link_libraries(p3plugin_add_tendencies PRIVATE p3plugin_global_vars p3plugin_types ComIn::ComIn)

target_link_libraries(p3plugin_main_routines PRIVATE p3_main p3plugin_utils p3plugin_global_vars p3plugin_types ComIn::ComIn)
target_link_libraries(${PROJECT_NAME} PRIVATE p3plugin_main_routines p3plugin_tracer_init p3plugin_utils p3plugin_global_vars ComIn::ComIn)
