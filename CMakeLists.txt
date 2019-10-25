cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

enable_language(Fortran)

set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90)

project(BlakeyFEM)

add_executable(test.exe ../src/test.f90)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)