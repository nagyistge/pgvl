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

ADD_EXECUTABLE( pgvl_cache_test
   CachePerformanceTest.cpp
)

#================Link======================

TARGET_LINK_LIBRARIES(pgvl_tests pgvl ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(pgvl_cache_test pgvl ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#================Tests=====================

ADD_TEST(
   NAME ImageTest
   COMMAND pgvl_tests --gtest_filter=ImageTest*
)

ADD_TEST(
   NAME ImageProcessingTest
   COMMAND pgvl_tests --gtest_filter=ImageProcessingTest*
)

#ADD_TEST(
#   NAME CachePerformanceTest
#   COMMAND valgrind --tool=callgrind --simulate-cache=yes --dump-instr=yes --collect-jumps=yes ${EXECUTABLE_OUTPUT_PATH}/pgvl_cache_test
#)
