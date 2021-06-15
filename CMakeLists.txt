cmake_minimum_required(VERSION 3.17)
#project(test)
project(lab11)

include_directories(.)

set(CMAKE_CXX_STANDARD 17)

#add_executable(${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}/${CMAKE_PROJECT_NAME}.cpp vector.h matrix.h)
#add_executable(${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}/kupa.cpp vector.h matrix.h)
#add_executable(${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}/kek.cpp vector.h matrix.h)

add_executable(${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}/${CMAKE_PROJECT_NAME}.cpp old/vector.h old/matrix.h old/vector.cpp old/matrix.cpp)
#add_executable(${CMAKE_PROJECT_NAME} lab5/lab5.cpp old2/vector_old.h old2/matrix_old.h)