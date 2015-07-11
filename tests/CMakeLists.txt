SET( PGVL_TEST_SRCS
   BitmapImageTest.cpp
   ImageProcessingTest.cpp
)

ADD_EXECUTABLE( pgvl_tests
   ${PGVL_TEST_SRCS}
)

# Fails to compile without pthread
FIND_PACKAGE(Threads REQUIRED)
TARGET_LINK_LIBRARIES(pgvl_tests pgvl ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

ADD_TEST(
   NAME BitmapImageTest
   COMMAND pgvl_tests --gtest_filter=BitmapImageTest*
)

ADD_TEST(
   NAME ImageProcessingTest
   COMMAND pgvl_tests --gtest_filter=ImageProcessingTest*
)