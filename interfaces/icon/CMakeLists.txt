cmake_minimum_required(VERSION 3.20)
project(p3_plugin LANGUAGES C Fortran)

set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
message(STATUS "Looking for FindNetCDF.cmake in: ${CMAKE_MODULE_PATH}")

# use this to fix for horeka module path bug while using lib/netcdf-fortran/4.6-intel-2025.1_llvm: 
set(NF_CONFIG_BIN "/software/all/lib/netcdf-fortran/4.6.2_intel_2025.1_llvm_serial/bin")

# or instead set to whatever lies in $NETCDF_FORTRAN_BIN_DIR or $NETCDF_FORTRAN_ROOT/bin (example):
#set(NF_CONFIG_BIN "/software/all/lib/netcdf-fortran/4.6_intel_2025.1_llvm_impi_2021.11/bin")

find_package(NetCDF REQUIRED COMPONENTS Fortran)

# compile the plugin with OpenMP for hybrid OpenMP-MPI parallelization:
#find_package(OpenMP REQUIRED COMPONENTS Fortran)

find_package(ComIn REQUIRED)

add_library(${PROJECT_NAME} MODULE p3_plugin.F90)
add_library(microphy_p3 SHARED microphy_p3.F90)


# compiler options without runtime optimizations and debugging options:
target_compile_options(${PROJECT_NAME} PRIVATE
  -g -traceback -debug -O0 -fma -fp-model precise -ftz -qno-opt-dynamic-align)
target_compile_options(microphy_p3 PRIVATE
  -g -traceback -debug -O0 -fma -fp-model precise -ftz -qno-opt-dynamic-align)

# compiler options with fast runtime optimizations without debugging:
#target_compile_options(${PROJECT_NAME} PRIVATE
#  -O2 -fma -fp-model fast -ftz -qno-opt-dynamic-align)
#target_compile_options(microphy_p3 PRIVATE
#  -O2 -fma -fp-model fast -ftz -qno-opt-dynamic-align)

target_link_libraries(${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran)

# compile the plugin with OpenMP for hybrid OpenMP-MPI parallelization:
#target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)

target_link_libraries(${PROJECT_NAME} PRIVATE microphy_p3 ComIn::ComIn)
