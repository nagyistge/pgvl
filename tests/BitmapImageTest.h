#ifndef BITMAPIMAGETEST_H
#define BITMAPIMAGETEST_H

#include <gtest/gtest.h>
#include "config.h"
#include "ppm.h"

class BitmapImageTest : public testing::Test {
public:
   BitmapImageTest();

   // From class Test
   virtual void SetUp();
   virtual void TearDown();

private:
};

// Test pgm loading
TEST_F(BitmapImageTest, loadsPgm) {
   int w = 0;
   int h = 0;
   unsigned char* data = pgmread(TEST_IMAGE_DIR "lena_gray.pgm", &w, &h);
   
   EXPECT_TRUE( data != 0 );
   EXPECT_EQ( w, 512 );
   EXPECT_EQ( h, 512 );
   
   delete[] data;
}

// Test ppm loading
TEST_F(BitmapImageTest, loadsPpm) {
   int w = 0;
   int h = 0;
   int maxval = 0;
   unsigned char* data = ppmread(TEST_IMAGE_DIR "lena.ppm", &w, &h, &maxval);
   
   EXPECT_TRUE( data != 0 );
   EXPECT_EQ( w, 512 );
   EXPECT_EQ( h, 512 );
   EXPECT_EQ( maxval, 255 );
   
   delete[] data;
}

#endif /*BITMAPIMAGETEST_H*/