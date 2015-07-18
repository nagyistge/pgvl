SET( PGVL_TEST_SRCS
   ImageTest.cpp
   ImageProcessingTest.cpp
)

# Fails to compile without pthread
FIND_PACKAGE(Threads REQUIRED)

#=============Executables==================

ADD_EXECUTABLE( pgvl_tests
   ${PGVL_TEST_SRCS}
)

IF( ${PERFORMANCE_TESTS} )
   ADD_EXECUTABLE( pgvl_cache_test
      CachePerformanceTest.cpp
   )
ENDIF()

#================Link======================

TARGET_LINK_LIBRARIES(pgvl_tests pgvl ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
IF( ${PERFORMANCE_TESTS} )
   TARGET_LINK_LIBRARIES(pgvl_cache_test pgvl ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

#================Tests=====================

ADD_TEST(
   NAME ImageTest
   COMMAND pgvl_tests --gtest_filter=ImageTest*
)

ADD_TEST(
   NAME ImageProcessingTest
   COMMAND pgvl_tests --gtest_filter=ImageProcessingTest*
)

IF( ${PERFORMANCE_TESTS} )
   ADD_TEST(
      NAME CachePerformanceTest
      COMMAND valgrind --tool=callgrind --simulate-cache=yes --dump-instr=yes --collect-jumps=yes ${EXECUTABLE_OUTPUT_PATH}/pgvl_cache_test
   )
ENDIF()
