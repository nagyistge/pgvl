#ifndef IMAGEPROCESSINGTEST_H
#define IMAGEPROCESSINGTEST_H

#include <Image.h>
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
   Image<uint8_t> img(rows, cols, chans);
   for( i = 0; i < rows; ++i )
      for( j = 0; j < cols*chans; ++j )
         img[i][j] = j + i*cols*chans;

   // Brute-force calculate expected integral image
   Image<uint8_t> expectedIntImg(img);
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
   Image<uint8_t> img(rows, cols, chans);
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

   Image<int8_t> img(2,3,1);
   img[0][0] = 1;
   img[0][1] = 2;
   img[0][2] = 3;
   img[1][0] = 9;
   img[1][1] = 6;
   img[1][2] = 5;

   Image<int8_t> kern(1,3,1);
   kern[0][0] = -1;
   kern[0][1] = 0;
   kern[0][2] = 1;

   Image<int8_t> out(2,3,1);

   filter(out, img, kern);
   EXPECT_EQ(out[0][1], 2);
   EXPECT_EQ(out[1][1], -4);
}

TEST_F(ImageProcessingTest, lowpassFilter) {
   Image<uint8_t> lena(TEST_IMAGE_DIR "lena_gray.pgm");

   Image<uint8_t> lpf(lena.rows(), lena.cols(), lena.channels());
   lowpassFilter(lpf, lena, 4);

   lpf.save("/tmp/lena_lpf");
}

// Tests the flow to rgb conversion, and also provides a flow key
TEST_F(ImageProcessingTest, opticalFlowToRgb) {
   int const radius = 32;
   int const rows = 2*radius+1;
   int const cols = rows;
   Image<float> flow(rows, cols, 2);
   Image<uint8_t> flowRgb(rows, cols, 3);

   for(int i = 0; i < rows; ++i) {
      for(int j = 0; j < cols; ++j) {
         flow[i][j*2+0] = static_cast<float>(j-radius)/radius;
         flow[i][j*2+1] = static_cast<float>(i-radius)/radius;
      }
   }

   opticalFlowToRgb(flowRgb, flow, 1.f);
   flowRgb.save("/tmp/flowkey");
}

TEST_F(ImageProcessingTest, hsOpticalFlow) {
   Image<uint8_t> frame1(TEST_IMAGE_DIR "rubic.0.pgm");
   Image<uint8_t> frame2(TEST_IMAGE_DIR "rubic.1.pgm");
   Image<float> fframe1;
   Image<float> fframe2;
   Image<float> flow(frame1.rows(), frame2.cols(), 2);
   Image<uint8_t> flowRgb(flow.rows(), flow.cols(), 3);

   fframe1.convertFrom(frame1);
   fframe2.convertFrom(frame2);
   hsOpticalFlow(flow, fframe1, fframe2);
   opticalFlowToRgb(flowRgb, flow);

   flowRgb.save("/tmp/flow");

   // Just for funsies, try out the SDL stuff ---------------------------------

   // Initialize SDL for video junk
   SDL_Init( SDL_INIT_VIDEO );

   // Create the window in which to display lena
   SDL_Window* window = SDL_CreateWindow(
      "Horn-Schunck Flow",
      SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
      frame1.cols(), frame1.rows(),
      SDL_WINDOW_SHOWN
   );
   // Create an accelerated renderer for the window
   SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
   // Convert Image to SDL_Texture for rendering
   SDL_Surface* imageSurface = toSurface(frame1);
   SDL_Texture* imageTexture = SDL_CreateTextureFromSurface(renderer, imageSurface);
   SDL_FreeSurface(imageSurface);

   SDL_RenderClear(renderer);
   SDL_RenderCopy(renderer, imageTexture, NULL, NULL);

   // Draw green lines indicating optical flow direction
   SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
   int const drawWidth = 8;
   float const flowScale = 10.f;
   for(int i = 0; i < frame1.rows(); i += drawWidth) {
      for(int j = 0; j < frame1.cols(); j += drawWidth) {
         SDL_RenderDrawLine(renderer, j, i, j+flowScale*flow[i][j*2+0]+0.5f, i+flowScale*flow[i][j*2+1]+0.5f);
      }
   }

   // Display for 10 seconds
   SDL_RenderPresent(renderer);
   SDL_Delay(10000);

   // Cleanup
   SDL_DestroyTexture(imageTexture);
   SDL_DestroyRenderer(renderer);
   SDL_DestroyWindow(window);
   SDL_Quit();
}

#endif /*IMAGEPROCESSINGTEST_H*/
