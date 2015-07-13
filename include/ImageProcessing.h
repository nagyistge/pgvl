/*
 * ImageProcessing.h is part of pgvl and is
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */

#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H

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
template<class T>
void filter(
   BitmapImage<T>& out,
   BitmapImage<T> const& img,
   BitmapImage<T> const& kernel,
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

#endif /*IMAGEPROCESSING_H*/