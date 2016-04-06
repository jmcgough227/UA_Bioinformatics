project(clustalw)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")
set(CMAKE_C_FLAGS "-lm")

add_executable(drawgram tree/draw.c tree/draw2.c tree/drawgram.c tree/phylip.c)
find_library(MATH m)
target_link_libraries(drawgram ${MATH})

add_executable(clustalw clustalw.cpp alignment_utilities.cpp parser.cpp)

target_link_libraries(clustalw quadmath)