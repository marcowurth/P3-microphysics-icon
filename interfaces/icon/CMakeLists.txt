cmake_minimum_required(VERSION 3.20)
project(p3_plugin LANGUAGES Fortran)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)

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

# compile the plugin with OpenMP for hybrid OpenMP-MPI parallelization:
#target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)

target_link_libraries(${PROJECT_NAME} PRIVATE microphy_p3 ComIn::ComIn)
