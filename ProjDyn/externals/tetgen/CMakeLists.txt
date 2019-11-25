file(GLOB SOURCES ./*.cpp)
file(GLOB HEADERS ./*.h)

if(UNIX)
    add_library(tetgen SHARED ${SOURCES} ${HEADERS})
elseif(WIN32)
    add_library(tetgen STATIC ${SOURCES} ${HEADERS})
endif()

set_target_properties(tetgen PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(tetgen PROPERTIES VERSION 1.0)
