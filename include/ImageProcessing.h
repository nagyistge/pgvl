/*
 * ImageProcessing.h is part of pgvl and is
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */

#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H

#include <cmath>
#include <BitmapImage.h>
#include <Point.h>

/*!
 * \defgroup ImageProcessing Image Processing
 * \brief All the basic image processing functionality
 */

/*!
 * \ingroup ImageProcessing
 * \brief Create an integral image
 *
 * For now, purposefully ignores overflowing
 *
 * \param[in,out] img the input image and output integral image
 */
template<class T>
void integrate(BitmapImage<T>& img) {
   int i,j;
   int const channels = img.channels();

   // Prefix scan all the rows
#pragma omp parallel for shared(img) private(i,j)
   for( i = 0; i < img.rows(); ++i ) {
      for( j = channels; j < img.rowWidth(); ++j )
         img[i][j] += img[i][j-channels];
   }

   // Prefix scan all the columns
#pragma omp parallel for shared(img) private(i,j)
   for( j = 0; j < img.rowWidth(); ++j ) {
      for( i = 1; i < img.rows(); ++i )
         img[i][j] += img[i-1][j];
   }
}

/*!
 * \ingroup ImageProcessing
 * \brief Create a squared integral image
 *
 * For now, purposefully ignores overflowing
 *
 * \param[in,out] img the input image and output squared integral image
 */
template<class T>
void integrateSquare(BitmapImage<T>& img) {
   int i,j;
   int const channels = img.channels();

   // Square the first pixel in each row
#pragma omp parallel for shared(img) private(i,j)
   for( i = 0; i < img.rows(); ++i )
      for( j = 0; j < channels; ++j )
         img[i][j] *= img[i][j];

   // Prefix scan all the rows
#pragma omp parallel for shared(img) private(i,j)
   for( i = 0; i < img.rows(); ++i ) {
      for( j = channels; j < img.rowWidth(); ++j )
         img[i][j] = (img[i][j]*img[i][j]) + img[i][j-channels];
   }

   // Prefix scan all the columns
   // NOTE: we do not have to square any pixels here, because
   // all pixels have been squared in the row scan above
#pragma omp parallel for shared(img) private(i,j)
   for( j = 0; j < img.rowWidth(); ++j ) {
      for( i = 1; i < img.rows(); ++i )
         img[i][j] += img[i-1][j];
   }
}

/*!
 * \ingroup ImageProcessing
 * \brief Filter an image with a kernel
 *
 * \note This does not do convolution, but correlation
 * \note For now, only does non-fft filtering
 *
 * \param out The output of the filtering
 * \param img Image to be filtered
 * \param kernel Kernel to apply to the \c img
 * \param anchor point within the kernel considered to be the center. If left
 *        at its default value, the anchor will be set to the center of the kernel.
 * \param delta value to add to the filtered value before storing in \c out
 */
template<class T, class U>
void filter(
   BitmapImage<T>& out,
   BitmapImage<T> const& img,
   BitmapImage<U> const& kernel,
   Point const& anchor = Point(-1,-1),
   float delta = 0.f
) {
   int const kcols = kernel.cols();
   int const krows = kernel.rows();

   int const rows = out.rows();
   int const cols = out.cols();
   int const channels = out.channels();

   int const anchorRow = (anchor==Point(-1,-1)) ? krows/2 : anchor.y;
   int const anchorCol = (anchor==Point(-1,-1)) ? kcols/2 : anchor.x;

   // Indices must always be valid:
   // img row: i + m - anchorRow
   // img col: j + n - anchorCol
   //
   // i + m - anchorRow >= 0:   i >= anchorRow
   // i + m - anchorRow < rows: i < rows + anchorRow - (krows-1)

   int i,j,k;
   int m,n;
#pragma omp parallel for shared(out,img,kernel) private(i,j,k,m,n)
   for( i = anchorRow; i < rows + anchorRow - krows + 1; ++i ) {
      for( j = anchorCol; j < cols + anchorCol - kcols + 1; ++j ) {
         for( k = 0; k < out.channels(); ++k ) {
            T& sum = out[i][j*channels + k];

            //sum = 0;
            // NOTE: is this the right place to apply the delta?
            sum = delta;
            for( m = 0; m < krows; ++m ) {
               for( n = 0; n < kcols; ++n ) {
                  sum += kernel[m][n*channels + k] * img[i+ m-anchorRow][(j+n-anchorCol)*channels + k];
               }
            }
         }
      }
   }
}

/*!
 * \ingroup ImageProcessing
 * \brief Specialization for uint8_t images
 */
template<>
void filter(
   BitmapImage<uint8_t>& out,
   BitmapImage<uint8_t> const& img,
   BitmapImage<uint8_t> const& kernel,
   Point const& anchor,
   float delta
) {
   int const kcols = kernel.cols();
   int const krows = kernel.rows();

   int const rows = out.rows();
   int const cols = out.cols();
   int const channels = out.channels();

   int const anchorRow = (anchor==Point(-1,-1)) ? krows/2 : anchor.y;
   int const anchorCol = (anchor==Point(-1,-1)) ? kcols/2 : anchor.x;

   // Indices must always be valid:
   // img row: i + m - anchorRow
   // img col: j + n - anchorCol
   //
   // i + m - anchorRow >= 0:   i >= anchorRow
   // i + m - anchorRow < rows: i < rows + anchorRow - (krows-1)

   int i,j,k;
   int m,n;
#pragma omp parallel for shared(out,img,kernel) private(i,j,k,m,n)
   for( i = anchorRow; i < rows + anchorRow - krows + 1; ++i ) {
      for( j = anchorCol; j < cols + anchorCol - kcols + 1; ++j ) {
         for( k = 0; k < out.channels(); ++k ) {
            // NOTE: is this the right place to apply the delta?
            int sum = delta;
            for( m = 0; m < krows; ++m ) {
               for( n = 0; n < kcols; ++n ) {
                  sum += static_cast<int>(kernel[m][n*channels + k]) * img[i+ m-anchorRow][(j+n-anchorCol)*channels + k];
               }
            }
            // Assume that 255 in the kernel corresponds to 1.0
            out[i][j*channels + k] = sum / 255;
         }
      }
   }
}

/*!
 * \ingroup ImageProcessing
 * \brief A 2D Gaussian function
 *
 * \tparam T the type of the x and y coordinates
 * \returns function float(T x, T y, T muX, T muY, T sigmaX, T sigmaY, T theta)
 *          where (x,y) are the coordinates, (muX,muY) are the means,
 *          (sigmaX,sigmaY) are the standard deviations, and theta is the angle
 */
template<class T>
std::function<float(T,T,T,T,T,T,float)> gauss() {
   return [](T x, T y, T muX, T muY, T sigmaX, T sigmaY, float theta) -> float {
      float const cosTheta = cosf(theta);
      float const sinTheta = sinf(theta);
      float const sin2Theta = sinf(2*theta);

      float const a = cosTheta*cosTheta/(2*sigmaX*sigmaX) + sinTheta*sinTheta/(2*sigmaY*sigmaY);
      float const b = -sin2Theta/(4*sigmaX*sigmaX) + sin2Theta/(4*sigmaY*sigmaY);
      float const c = sinTheta*sinTheta/(2*sigmaX*sigmaX) + cosTheta*cosTheta/(2*sigmaY*sigmaY);

      return expf( -(a*(x-muX)*(x-muX) + 2*b*(x-muX)*(y-muY) + c*(y-muY)*(y-muY)) );
   };
}

/*!
 * \ingroup ImageProcessing
 * \brief Image lowpass filtering
 *
 * \param[out] out output image
 * \param[in] img input image
 * \param[in] radius spatial radius in pixels of the lowpass filter
 */
template<class T>
void lowpassFilter(
   BitmapImage<T>& out,
   BitmapImage<T> const& img,
   int radius
) {
   auto kFunc = gauss<int>();
   int const kSize = radius % 2 == 0 ? 3*radius+1 : 3*radius;
   int const kCenter = kSize/2;
   int const kStd = radius/2;

   int const chans = img.channels();
   BitmapImage<float> kernelX(1,kSize,img.channels());
   BitmapImage<float> kernelY(kSize,1,img.channels());
   float sum = 0.f;
   float val;
   for(int i = 0; i < kSize; ++i) {
      val = kFunc(i, 0, kCenter, 0, kStd, 1, 0.f);
      sum += val;
      for(int k = 0; k < chans; ++k)
         kernelX[0][i*chans+k] = kernelY[i][k] = val;
   }
   for(int i = 0; i < kSize; ++i) {
      for(int k = 0; k < chans; ++k) {
         kernelX[0][i*chans+k] /= sum;
         kernelY[i][k] /= sum;
      }
   }

   BitmapImage<T> tmp(img.rows(), img.cols(), img.channels());
   filter(tmp, img, kernelX);
   filter(out, tmp, kernelY);
}

/*!
 * \ingroup ImageProcessing
 * \brief Specialization for uint8_t images
 */
template<>
void lowpassFilter(
   BitmapImage<uint8_t>& out,
   BitmapImage<uint8_t> const& img,
   int radius
) {
   auto kFunc = gauss<int>();
   int const kSize = radius % 2 == 0 ? 3*radius+1 : 3*radius;
   int const kCenter = kSize/2;
   int const kStd = radius/2;

   int const chans = img.channels();
   BitmapImage<uint8_t> kernelX(1,kSize,img.channels());
   BitmapImage<uint8_t> kernelY(kSize,1,img.channels());
   float sum = 0;
   int val;
   for(int i = 0; i < kSize; ++i) {
      val = 255.f * kFunc(i, 0, kCenter, 0, kStd, 1, 0.f);
      sum += val;
      for(int k = 0; k < chans; ++k)
         kernelX[0][i*chans+k] = kernelY[i][k] = val;
   }
   sum /= 255.f;
   for(int i = 0; i < kSize; ++i) {
      for(int k = 0; k < chans; ++k) {
         kernelX[0][i*chans+k] /= sum;
         kernelY[i][k] /= sum;
      }
   }

   BitmapImage<uint8_t> tmp(img.rows(), img.cols(), img.channels());
   filter(tmp, img, kernelX);
   filter(out, tmp, kernelY);
}

#endif /*IMAGEPROCESSING_H*/