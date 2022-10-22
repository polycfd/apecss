cmake_minimum_required(VERSION 3.12)

project (gastemperature_apecss)
set(CMAKE_C_COMPILER /usr/bin/gcc)

include_directories($ENV{APECSS_DIR}/include)
FILE (GLOB_RECURSE myfiles ABSOLUTE ../src/*.c)

set (mylibs m apecss)
link_directories($ENV{USRLIB_DIR} $ENV{APECSS_DIR}/lib)

foreach(arg ${myincludes})
  IF (arg MATCHES "-I")
    STRING(REGEX REPLACE "-I" "" myinc ${arg})
    message("Additional include: ${myinc}")
    include_directories(${myinc})
  ENDIF(arg MATCHES "-I")
endforeach(arg ${myincludes})

foreach(arg ${mylibs})
  STRING(REGEX REPLACE "lib" "" myl1 ${arg})
  STRING(REGEX REPLACE ".a$" "" myl2 ${myl1})
  set(mylibs ${myl2} ${mylibs} ${myl2})
endforeach(arg ${mylibs})

if(NOT CMAKE_BUILD_TYPE)
   message(STATUS "Optimization: No optimization specified.")
   set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
   message(STATUS "Optimization: Debug")
   set(CMAKE_C_FLAGS_DEBUG "-Wall -g")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
   message(STATUS "Optimization: Release")
   add_definitions(-DNDEBUG)
   set(CMAKE_C_FLAGS_RELEASE "-Wall -Werror -O3")
endif()

add_executable(gastemperature_apecss ${myfiles})
target_link_libraries(gastemperature_apecss ${mylibs})