set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)

PROJECT(coortrans)
include_directories(include)
AUX_SOURCE_DIRECTORY(src DIR_SRCS)
SET(TEST_MATH ${DIR_SRCS})
add_executable(a ${TEST_MATH})