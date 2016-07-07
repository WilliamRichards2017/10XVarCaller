INCLUDE_DIRECTORIES(
  ${BAMTOOLS_INCLUDE}
  ${ZLIB_INCLUDE}
  )

set(10X_TOOLS_SOURCES
  main.cpp
)

add_executable(10x 
${10X_TOOLS_SOURCES}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


add_dependencies(10x ${ZLIB_PROJECT} ${BAMTOOLS_PROJECT})

target_link_libraries(10x
    ${BAMTOOLS_LIB}
    ${BAMTOOLS_UTIL_LIB}
    ${ZLIB_LIB}
 )
