#ifndef IMAGEPROCESSINGTEST_H
#define IMAGEPROCESSINGTEST_H

#include <BitmapImage.h>
#include <ImageProcessing.h>
#include <gtest/gtest.h>
#include <stdio.h>

class ImageProcessingTest : public testing::Test {
public:
   ImageProcessingTest();
   virtual void SetUp();
   virtual void TearDown();
private:
};

TEST_F(ImageProcessingTest, integrate) {
   int i,j;
   int const rows = 2;
   int const cols = 3;
   int const chans = 3;

   // img:
   // [ 0  1  2 |  3  4  5 |  6  7  8
   //   9 10 11 | 12 13 14 | 15 16 17
   //  18 19 20 | 21 22 23 | 24 25 26]
   BitmapImage img(rows, cols, chans);
   for( i = 0; i < rows; ++i )
      for( j = 0; j < cols*chans; ++j )
         img[i][j] = j + i*cols*chans;

   // Brute-force calculate expected integral image
   BitmapImage expectedIntImg(img);
   for( j = chans; j < cols*chans; ++j )
      expectedIntImg[0][j] += expectedIntImg[0][j-chans];
   for( i = 1; i < rows; ++i ) {
      for( j = 0; j < chans; ++j )
         expectedIntImg[i][j] += expectedIntImg[i-1][j];
   }
   for( i = 1; i < rows; ++i )
      for( j = chans; j < cols*chans; ++j )
         expectedIntImg[i][j] += expectedIntImg[i-1][j] + expectedIntImg[i][j-chans] - expectedIntImg[i-1][j-chans];

   // Test
   integrate(img);
   for( i = 0; i < rows; ++i )
      for( j = 0; j < cols*chans; ++j )
         EXPECT_EQ( img[i][j], expectedIntImg[i][j] );
}

#endif /*IMAGEPROCESSINGTEST_H*/