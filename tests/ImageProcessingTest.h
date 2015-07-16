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
   BitmapImage<uint8_t> img(rows, cols, chans);
   for( i = 0; i < rows; ++i )
      for( j = 0; j < cols*chans; ++j )
         img[i][j] = j + i*cols*chans;

   // Brute-force calculate expected integral image
   BitmapImage<uint8_t> expectedIntImg(img);
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

TEST_F(ImageProcessingTest, integrateSquare) {
   int i,j;
   int const rows = 2;
   int const cols = 3;
   int const chans = 1;

   // img:
   // [ 0  1  2
   //   3  4  5 ]
   BitmapImage<uint8_t> img(rows, cols, chans);
   for( i = 0; i < rows; ++i )
      for( j = 0; j < cols*chans; ++j )
         img[i][j] = j + i*cols*chans;

   // img^2:
   // [ 0  1  4
   //   9 16 25 ]

   // Row prefix
   // [ 0  1  5
   //   9 25 50 ]

   // Col prefix
   // [ 0  1  5
   //   9 26 55 ]
   integrateSquare(img);
   EXPECT_EQ( 0 + 1 + 4 + 9 + 16 + 25, img[1][2] );
}

TEST_F(ImageProcessingTest, filter) {

   // img:
   // [ 1 2 3
   //   9 6 5 ]

   // kern:
   // [ -1 0 1 ]

   // out:
   // [ 0  2 0
   //   0 -4 0 ]

   BitmapImage<int8_t> img(2,3,1);
   img[0][0] = 1;
   img[0][1] = 2;
   img[0][2] = 3;
   img[1][0] = 9;
   img[1][1] = 6;
   img[1][2] = 5;

   BitmapImage<int8_t> kern(1,3,1);
   kern[0][0] = -1;
   kern[0][1] = 0;
   kern[0][2] = 1;

   BitmapImage<int8_t> out(2,3,1);

   filter(out, img, kern);
   EXPECT_EQ(out[0][1], 2);
   EXPECT_EQ(out[1][1], -4);
}

TEST_F(ImageProcessingTest, lowpassFilter) {
   BitmapImage<uint8_t> lena(TEST_IMAGE_DIR "lena_gray.pgm");

   BitmapImage<uint8_t> lpf(lena.rows(), lena.cols(), lena.channels());
   lowpassFilter(lpf, lena, 4);

   lpf.save("/tmp/lena_lpf");
}

TEST_F(ImageProcessingTest, hsOpticalFlow) {
   BitmapImage<uint8_t> frame1(TEST_IMAGE_DIR "office.0.ppm");
   BitmapImage<uint8_t> frame2(TEST_IMAGE_DIR "office.1.ppm");
   BitmapImage<float> fframe1;
   BitmapImage<float> fframe2;
   BitmapImage<float> flow(frame1.rows(), frame2.cols(), 2);
   BitmapImage<uint8_t> flowRgb(flow.rows(), flow.cols(), 3);

   fframe1.convertFrom(frame1);
   fframe2.convertFrom(frame2);
   hsOpticalFlow(flow, fframe1, fframe2);
   opticalFlowToRgb(flowRgb, flow);

   flowRgb.save("/tmp/flow");
}

#endif /*IMAGEPROCESSINGTEST_H*/