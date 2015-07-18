/*
 * ImageProcessing.h is part of pgvl and is
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */

#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H

#include <cmath>
#include <Image.h>
#include <Point.h>
#include <Eigen/Dense>

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
void integrate(Image<T>& img) {
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
void integrateSquare(Image<T>& img) {
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
   Image<T>& out,
   Image<T> const& img,
   Image<U> const& kernel,
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
         for( k = 0; k < channels; ++k ) {
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
   Image<uint8_t>& out,
   Image<uint8_t> const& img,
   Image<uint8_t> const& kernel,
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
         for( k = 0; k < channels; ++k ) {
            // NOTE: is this the right place to apply the delta?
            int sum = delta;
            for( m = 0; m < krows; ++m ) {
               for( n = 0; n < kcols; ++n ) {
                  int a = kernel[m][n*channels + k];
                  // TODO: this line causes a lot of L1 misses
                  int b = img[i+m-anchorRow][(j+n-anchorCol)*channels + k];
                  sum += a*b;
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
   Image<T>& out,
   Image<T> const& img,
   int radius
) {
   auto kFunc = gauss<int>();
   int const kSize = radius % 2 == 0 ? 3*radius+1 : 3*radius;
   int const kCenter = kSize/2;
   int const kStd = radius/2;

   int const chans = img.channels();
   Image<float> kernelX(1,kSize,img.channels());
   Image<float> kernelY(kSize,1,img.channels());
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

   Image<T> tmp(img.rows(), img.cols(), img.channels());
   filter(tmp, img, kernelX);
   filter(out, tmp, kernelY);
}

/*!
 * \ingroup ImageProcessing
 * \brief Specialization for uint8_t images
 */
template<>
void lowpassFilter(
   Image<uint8_t>& out,
   Image<uint8_t> const& img,
   int radius
) {
   auto kFunc = gauss<int>();
   int const kSize = radius % 2 == 0 ? 3*radius+1 : 3*radius;
   int const kCenter = kSize/2;
   int const kStd = radius/2;

   int const chans = img.channels();
   Image<uint8_t> kernelX(1,kSize,img.channels());
   Image<uint8_t> kernelY(kSize,1,img.channels());
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

   Image<uint8_t> tmp(img.rows(), img.cols(), img.channels());
   filter(tmp, img, kernelX);
   filter(out, tmp, kernelY);
}

/*!
 * \ingroup ImageProcessing
 * \brief Get spatial gradients
 */
template<class T>
void gradient(
   Image<float>& outDx,
   Image<float>& outDy,
   Image<T> const& img
) {

   int const chans = img.channels();
   Image<float> dx(1,3,img.channels());
   Image<float> dy(3,1,img.channels());

   dx[0][0*chans + 0] = -0.5f;
   dx[0][2*chans + 0] = 0.5f;

   dy[0][0] = -0.5f;
   dy[2][0] = 0.5f;

   for(int i = 1; i < chans; ++i) {
      dx[0][0*chans + i] = dx[0][0*chans + 0];
      dx[0][2*chans + i] = dx[0][2*chans + 0];
      dy[0][chans] = dy[0][0];
      dy[2][chans] = dy[2][0];
   }

   filter(outDx, img, dx);
   filter(outDy, img, dy);
}

/*!
 * \ingroup ImageProcessing
 * \brief Convert dense optical flow to an rgb image for display
 *
 * \param[out] rgb output rgb image
 * \param[in] flow flow image whose 2 channels are [vx, vy]
 * \param[in] maxFlow denominator used to scale flow values for display
 */
void opticalFlowToRgb(
   Image<uint8_t>& rgb,
   Image<float>& flow,
   float const maxFlow = 2.f
){
   int const rows = rgb.rows();
   int const cols = rgb.cols();

   Image<float> hsv(rows, cols, 3);
   int i,j;
   for( i = 0; i < rows; ++i ) {
      for( j = 0; j < cols; ++j ) {
         float dx = flow[i][j*2+0];
         float dy = flow[i][j*2+1];
         float angle = 180.f/M_PI * atan2(dy, dx);
         float mag = sqrtf(dx*dx + dy*dy) / maxFlow;

         if( angle < 0.f )
            angle += 360.f;

         hsv[i][j*3+0] = angle;
         hsv[i][j*3+1] = std::min(1.f, mag);
         hsv[i][j*3+2] = 1.f;
      }
   }

   hsv2rgb(hsv);
   rgb.convertFrom<float>(hsv, [](float x) -> uint8_t { return static_cast<uint8_t>(255.f*x); });
}

/*!
 * \ingroup ImageProcessing
 * \brief Horn-Schunck optical flow
 *
 * \param[out] flow output optical flow image, whose 2 channels are [vx, vy]
 * \param[in] img0 reference image frame
 * \param[in] img1 image frame coming temporally after \c img0
 */
void hsOpticalFlow(
   Image<float>& flow,
   Image<float> const& img0,
   Image<float> const& img1
) {
   // Patch radius in pixels
   int const radius = 3;
   int const rows = img0.rows();
   int const cols = img0.cols();
   int const chans = img0.channels();
   // Regularization parameter (bias flow towards 0)
   float const gamma = 1e-2 * (2*radius+1)*(2*radius+1)*chans;
   int i,j,k;
   int m,n;
   Image<float> x0(rows, cols, chans);
   Image<float> x1(rows, cols, chans);
   Image<float> dx(rows, cols, chans);
   Image<float> dy(rows, cols, chans);
   Image<float> dt(rows, cols, chans);

   // LPF for noise
   lowpassFilter(x0, img0, 2);
   lowpassFilter(x1, img1, 2);

   // Get dx,dy,dt
   gradient(dx, dy, x0);
   for(i = 0; i < rows; ++i)
      for(j = 0; j < cols; ++j)
         for(k = 0; k < chans; ++k)
            dt[i][j*chans+k] = x0[i][j*chans+k] - x1[i][j*chans+k];

   // | Ix[0]  Iy[0] | [dx;dy] = | It[0] |
   // | Ix[1]  Iy[1] |           | It[1] |
   //   ...
   // | Ix[n]  Iy[n] |           | It[n] |
   //
   // A x = b
   // x = (A^TA)^-1 A^Tb
   // A^TA = [ \sum_{i=1}^n Ix[i]^2 & \sum{i=1}^n Ix[i]Iy[i] \\ \sum{i=1}^n Ix[i]Iy[i] & \sum{i=1}^n Iy[i]^2 ]
   // A^Tb = [ \sum_{i=1}^n Ix[i]It[i] \\ \sum_{i=1}^n Iy[i]It[i] ]

   Eigen::Matrix2f A;
   Eigen::Vector2f b;
   Eigen::Vector2f x;
   for(i = radius; i < rows - radius; ++i) {
      for(j = radius; j < cols - radius; ++j) {
         A.setZero();
         b.setZero();

         // NOTE: optimize this aggregation later with integral images
         for(m = i-radius; m < i+radius; ++m) {
            for(n = j-radius; n < j+radius; ++n) {
               for(k = 0; k < chans; ++k) {
                  A(0,0) += dx[m][n*chans+k]*dx[m][n*chans+k];
                  A(0,1) += dx[m][n*chans+k]*dy[m][n*chans+k];
                  A(1,1) += dy[m][n*chans+k]*dy[m][n*chans+k];

                  b(0) += dx[m][n*chans+k]*dt[m][n*chans+k];
                  b(1) += dy[m][n*chans+k]*dt[m][n*chans+k];
               }
            }
         }
         // Make it symmetric
         A(1,0) = A(0,1);
         // Apply regularization
         A(0,0) += gamma;
         A(1,1) += gamma;

         x = A.ldlt().solve(b);
         flow[i][j*2+0] = x(0);
         flow[i][j*2+1] = x(1);
      }
   }
}

#endif /*IMAGEPROCESSING_H*/
