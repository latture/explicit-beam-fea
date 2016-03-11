set(EXPLICIT_FEA_SOURCES
        beam_element.cpp
        euler_bernoulli_beam_element.cpp
        explicit_system.cpp
        fe_utils.cpp
        mesh.cpp
        setup.cpp
        timoshenko_beam_element.cpp
        explicit_system_manager.cpp)

INCLUDE_DIRECTORIES(/usr/include)
INCLUDE_DIRECTORIES(/opt/intel/include)
LINK_DIRECTORIES(/usr/lib)
LINK_DIRECTORIES(/opt/intel/lib)

add_library(explicit_fea ${EXPLICIT_FEA_SOURCES} ${BOOST_SOURCES})
add_executable(explicit_fea_cmd main.cpp)
target_link_libraries(explicit_fea_cmd explicit_fea pthread m dl)

if(EXPLICIT_FEA_BUILD_UNIT_TESTS)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # force this option to ON so that Google Test will use /MD instead of /MT
    # /MD is now the default for Visual Studio, so it should be our default, too
    option(gtest_force_shared_crt
           "Use shared (DLL) run-time lib even when Google Test is built as static lib."
           ON)
    endif()

    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

    MACRO (EXPLICIT_FEA_TEST NAME)
        add_executable(${NAME}_test ${NAME}_test.cpp)
        target_link_libraries(${NAME}_test explicit_fea gtest gtest_main)
        add_test(name ${NAME}_test
                 command ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME}_test
                 --test_srcdir
                 ${CMAKE_SOURCE_DIR}/data)
    ENDMACRO (EXPLICIT_FEA_TEST)

    EXPLICIT_FEA_TEST(euler_bernoulli_beam_element)
    EXPLICIT_FEA_TEST(explicit_system)
    EXPLICIT_FEA_TEST(explicit_system_manager)
    EXPLICIT_FEA_TEST(mesh)
    EXPLICIT_FEA_TEST(setup)
    EXPLICIT_FEA_TEST(timoshenko_beam_element)

endif(EXPLICIT_FEA_BUILD_UNIT_TESTS)