/*
 * ImageProcessing.h is part of pgvl and is
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */

#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H

#include <BitmapImage.h>

/*!
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

#endif /*IMAGEPROCESSING_H*/