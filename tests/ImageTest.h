#ifndef BITMAPIMAGETEST_H
#define BITMAPIMAGETEST_H

#include <gtest/gtest.h>
#include "config.h"
#include <ppm.h>
#include <Image.h>

class ImageTest : public testing::Test {
public:
   ImageTest();

   // From class Test
   virtual void SetUp();
   virtual void TearDown();

private:
};

// Test pgm loading
TEST_F(ImageTest, loadsPgm) {
   int w = 0;
   int h = 0;
   unsigned char* data = pgmread(TEST_IMAGE_DIR "lena_gray.pgm", &w, &h);

   EXPECT_TRUE( data != 0 );
   EXPECT_EQ( w, 512 );
   EXPECT_EQ( h, 512 );

   delete[] data;
}

// Test ppm loading
TEST_F(ImageTest, loadsPpm) {
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

// Test Image loading
TEST_F(ImageTest, loadsImage) {
   Image<uint8_t> ppm(TEST_IMAGE_DIR "lena.ppm");

   EXPECT_EQ( ppm.rows(), 512 );
   EXPECT_EQ( ppm.cols(), 512 );

   Image<uint8_t> pgm(TEST_IMAGE_DIR "lena_gray.pgm");

   EXPECT_EQ( pgm.rows(), 512 );
   EXPECT_EQ( pgm.cols(), 512 );
}

// Ensure saving then loading results in the same data
TEST_F(ImageTest, loadStore) {
   int i,j,k;
   Image<uint8_t> lena(TEST_IMAGE_DIR "lena_gray.pgm");
   lena.save("/tmp/pgvl-lena");
   Image<uint8_t> lena2("/tmp/pgvl-lena.pgm");

   EXPECT_EQ( lena.rows(), lena2.rows() );
   EXPECT_EQ( lena.cols(), lena2.cols() );
   EXPECT_EQ( lena.channels(), lena2.channels() );

   bool different = false;
   for( i = 0; i < lena.rows() && !different; ++i ) {
      for( j = 0; j < lena.cols(); ++j ) {
         if( lena[i][j] != lena2[i][j] ) {
            different = true;
            break;
         }
      }
   }

   if( different )
      std::cerr << "(" << i << ", " << j << ")" << std::endl;
   EXPECT_EQ( false, different );

   Image<uint8_t> lenaColor(TEST_IMAGE_DIR "lena.ppm");
   lenaColor.save("/tmp/pgvl-lena");
   Image<uint8_t> lena3("/tmp/pgvl-lena.ppm");

   different = false;
   for( i = 0; i < lenaColor.rows() && !different; ++i ) {
      for( j = 0; j < lenaColor.cols(); ++j ) {
         for( k = 0; k < lenaColor.channels(); ++k ) {
            if( lenaColor[i][j*lenaColor.channels()+k] != lena3[i][j*lenaColor.channels()+k] ) {
               different = true;
               break;
            }
         }
      }
   }

   if( different )
      std::cerr << "(" << i << ", " << j << ", " << k << ")" << std::endl;
   EXPECT_EQ( false, different );
}

TEST_F(ImageTest, patch) {
   Image<uint8_t> lena(TEST_IMAGE_DIR "lena_gray.pgm");
   Image<uint8_t> patch;

   lena.patch(patch, 1, 5, 1, 5);

   EXPECT_EQ( patch.rows(), 5 );
   EXPECT_EQ( patch.cols(), 5 );
   EXPECT_EQ( patch.channels(), lena.channels() );

   bool different = false;
   for( int i = 0; i < patch.rows() && !different; ++i ) {
      for( int j = 0; j < patch.cols(); ++j ) {
         if( patch[i][j] != lena[1+i][1+j] ) {
            different = true;
            break;
         }
      }
   }

   EXPECT_EQ( different, false );
}

#endif /*BITMAPIMAGETEST_H*/
