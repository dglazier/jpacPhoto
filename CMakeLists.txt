cmake_minimum_required (VERSION 2.9)


project (jpacPhoto)

set(CMAKE_CXX_FLAGS "-fPIC -O3") 

#Compare the new contents with the existing file, if it exists and is the
# same we don't want to trigger a make by changing its timestamp.
function(update_file path content)
   set(old_content "")
   if(EXISTS "${path}")
       file(READ "${path}" old_content)
   endif()
   if(NOT old_content STREQUAL content)
       file(WRITE "${path}" "${content}")
   endif()
endfunction(update_file)

# Creates a file called CMakeDeps.cmake next to your CMakeLists.txt with
# the list of dependencies in it - this file should be treated as part of
# CMakeLists.txt (source controlled, etc.).
function(update_deps_file deps)
    set(deps_file "CMakeDeps.cmake")
    # Normalize the list so it's the same on every machine
    list(REMOVE_DUPLICATES deps)
    foreach(dep IN LISTS deps)
        file(RELATIVE_PATH rel_dep ${CMAKE_CURRENT_SOURCE_DIR} ${dep})
        list(APPEND rel_deps ${rel_dep})
    endforeach(dep)
    list(SORT rel_deps)
    # Update the deps file
    set(content "# generated by make process\nset(sources ${rel_deps})\n")
    update_file(${deps_file} "${content}")
    # Include the file so it's tracked as a generation dependency we don't
    # need the content.
    include(${deps_file})
endfunction(update_deps_file)

# Look for and connect ROOT
find_package(ROOT REQUIRED MathMore)
if (ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIRS})
    link_directories(${ROOT_LIBRARY_DIRS})
else()
    message(SEND_ERROR "Root not found!")
endif()

# BUILD LIBRARY FROM LOCAL FiLES
include_directories("include")
include_directories("src")
file(GLOB_RECURSE INC "include/*.hpp")
file(GLOB_RECURSE SRC "src/*.cpp")
add_library( jpacPhoto SHARED ${INC} ${SRC} )

# BUILD THE PLOTTING LIBRARY AS WELL
include_directories("jpacStyle/include")
include_directories("jpacStyle/src")
file(GLOB_RECURSE PLOTINC "jpacStyle/include/*.hpp")
file(GLOB_RECURSE PLOTSRC "jpacStyle/src/*.cpp")
add_library( jpacStyle SHARED ${PLOTINC} ${PLOTSRC} )

# INSTALLATION SETTINGS
set( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

install(TARGETS jpacPhoto jpacStyle
  LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}" )

# GRAB EXECUTABLES and add to the makelist
include_directories("executables")
file(GLOB_RECURSE EXE_FILES "executables/*.cpp")
update_deps_file("${EXE_FILES}")
foreach( exefile ${EXE_FILES} )
    get_filename_component( exename ${exefile} NAME_WE)
    add_executable( ${exename} ${exefile} ) # MAKE SURE MAKE INSTALL DOESNT MAKE ALL THE EXECUTABLES
    target_link_libraries( ${exename} jpacPhoto)
    target_link_libraries( ${exename} jpacStyle)
    target_link_libraries( ${exename} ${ROOT_LIBRARIES})
endforeach( exefile ${EXE_FILES} )