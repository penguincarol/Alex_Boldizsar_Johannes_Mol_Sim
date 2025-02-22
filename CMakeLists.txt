#############################################################################################################
#   General                                                                                                 #
#############################################################################################################
cmake_minimum_required(VERSION 3.10)

# define project name, version
project(PSEMolDyn_GroupF VERSION 0.0.1)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# collect all cpp files
file(GLOB_RECURSE MY_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")

#############################################################################################################
#   Create Executable                                                                                       #
#############################################################################################################
add_subdirectory(libs/libspdlog)
#find_package(spdlog REQUIRED PATHS ${CMAKE_CURRENT_SOURCE_DIR}/libs/libspdlog)
add_executable(MolSim ${MY_SRC})
target_compile_features(MolSim PRIVATE cxx_std_20)
target_include_directories(MolSim
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libeigen/Eigen
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libspdlog/include/spdlog
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

find_package(PkgConfig REQUIRED)
pkg_check_modules(XERCESC REQUIRED IMPORTED_TARGET xerces-c)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MolSim
            PUBLIC PkgConfig::XERCESC
            PRIVATE spdlog::spdlog
            PUBLIC OpenMP::OpenMP_CXX
            )
else()
    target_link_libraries(MolSim
            PUBLIC PkgConfig::XERCESC
            PRIVATE spdlog::spdlog
            )
endif()

#this code got inflated really quickly..
# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
if(three_dim_tasks)
    if(round_robin_distr)
        target_compile_options(MolSim
                PRIVATE
                $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas -O3 -march=native -DTHREE_DIMENSIONAL_TASKS -DTASK_ROUND_ROBIN>
                # disable some spam warnings for icpc...
                $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
                )
    else()
        target_compile_options(MolSim
                PRIVATE
                $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas -O3 -march=native -DTHREE_DIMENSIONAL_TASKS>
                # disable some spam warnings for icpc...
                $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
                )
    endif()
elseif(one_dim_tasks)
    target_compile_options(MolSim
            PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas -O3 -march=native -DONE_DIMENSIONAL_TASKS>
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
            )
else()
    if(round_robin_distr)
        target_compile_options(MolSim
                PRIVATE
                $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas -O3 -march=native -DTASK_ROUND_ROBIN>
                # disable some spam warnings for icpc...
                $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
                )
    else()
        if(task_oriented_2d)
            target_compile_options(MolSim
                    PRIVATE
                    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas -O3 -march=native -DTASK_ORIENTED_2D>
                    # disable some spam warnings for icpc...
                    $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
                    )
        else()
        target_compile_options(MolSim
                PRIVATE
                $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas -O3 -march=native>
                # disable some spam warnings for icpc...
                $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
                )
        endif()
    endif()
endif()

#Benchmark option
#this is almost the identical code to the target MolSim and i am sure that there is an extremely pretty way to do exactely that in ~4 lines instead
add_executable(ProfileMolSim ${MY_SRC})
set_target_properties(ProfileMolSim PROPERTIES EXCLUDE_FROM_ALL True)
target_compile_features(ProfileMolSim PRIVATE cxx_std_20)
target_include_directories(ProfileMolSim
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libeigen/Eigen
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libspdlog/include/spdlog
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

find_package(PkgConfig REQUIRED)
pkg_check_modules(XERCESC REQUIRED IMPORTED_TARGET xerces-c)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ProfileMolSim
            PUBLIC PkgConfig::XERCESC
            PRIVATE spdlog::spdlog
            PUBLIC OpenMP::OpenMP_CXX
            )
else()
    target_link_libraries(ProfileMolSim
            PUBLIC PkgConfig::XERCESC
            PRIVATE spdlog::spdlog
            )
endif()
# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(ProfileMolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas -pg>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>

        )
target_link_libraries(ProfileMolSim PRIVATE -pg)

#copy format file
add_custom_command(TARGET MolSim PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different  "${CMAKE_CURRENT_SOURCE_DIR}/src/io/input/sim_input/XMLFormat.xsd" "${CMAKE_CURRENT_SOURCE_DIR}/resources/XMLFormat.xsd")
add_custom_command(TARGET MolSim PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different  "${CMAKE_CURRENT_SOURCE_DIR}/src/io/input/sim_input/XMLFormat.xsd" "${CMAKE_CURRENT_SOURCE_DIR}/resources/tests/XMLFormat.xsd")
add_custom_command(TARGET MolSim PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different  "${CMAKE_CURRENT_SOURCE_DIR}/src/io/input/sim_input/XMLFormat.xsd" "${CMAKE_CURRENT_SOURCE_DIR}/input/XMLFormat.xsd")
add_custom_command(TARGET MolSim POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/resources" "${CMAKE_BINARY_DIR}")
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# create sanitization target
add_executable(asan-MolSim ${MY_SRC})
set_target_properties(asan-MolSim PROPERTIES EXCLUDE_FROM_ALL True)
target_compile_features(asan-MolSim PRIVATE cxx_std_20)
target_include_directories(asan-MolSim
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libeigen/Eigen
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libspdlog/include/spdlog
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
        )
target_link_libraries(asan-MolSim PUBLIC PkgConfig::XERCESC PRIVATE spdlog::spdlog PUBLIC OpenMP::OpenMP_CXX)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(asan-MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall -Wno-unknown-pragmas>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
        PRIVATE -fsanitize=address
        )
target_link_options(asan-MolSim PRIVATE -fsanitize=address)


#############################################################################################################
#   GOOGLETEST                                                                                              #
#############################################################################################################
enable_testing()
include(FetchContent)
FetchContent_Declare(
        googletest
        # Specify the commit you depend on and update it regularly.
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(tests EXCLUDE_FROM_ALL)

#############################################################################################################
#   DOXYGEN                                                                                                 #
#############################################################################################################
include(doxygen)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Source: https://vicrucann.github.io/tutorials/quick-cmake-doxygen/
option(BUILD_DOC "Build documentation" ON)
find_package(Doxygen)

if(BUILD_DOC)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif(BUILD_DOC)
