/* ppm.h is part of pgvl and is 
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */
#ifndef PGVLPPM_H
#define PGVLPPM_H

/*!
 * \brief Reads a binary or ascii PGM (grayscale image) file.
 * 
 * \param filename The pgm file to read
 * \param w Output width
 * \param h Output height
 * \returns a pointer to the image in row-major format. The calling
 *          function has responsibility to delete[] the memory. A NULL is returned
 *          in the case of failure.
 */
unsigned char* pgmread(const char* filename, int* w, int* h);

/*!
 * \brief Read a normalized grayscale floating-point image.
 * 
 * \sa pgmread()
 */
float* pgmread_float(const char* filename, int* w, int* h );

unsigned char* rgbToArgb(unsigned char* rgb, int w, int h);

/*!
 * \brief Reads a binary or ascii PPM (rgb image) file.
 * 
 * \param filename The ppm file to read
 * \param w Output width
 * \param h Output height
 * \param maxval Output maximum value.
 * \returns a pointer to the image in row-major RGB format. The calling
 *          function has responsibility to delete[] the memory. A NULL is returned
 *          in the case of failure.
 */
unsigned char* ppmread(const char* filename, int* w, int* h, int* maxval);

/*!
 * \brief Read a normalized floating-point image.
 * 
 * \sa ppmread()
 */
float* ppmread_float(const char* filename, int* w, int* h );

/*!
 * \brief Write a PGM image.
 * 
 * \param filename The file to write to.
 * \param w Image width
 * \param h Image height
 * \param pitch number of bytes in one row of data
 * \param data Row-major image data
 * \param comment_string Comments (NULL if none)
 * \param binsave true for binary writing, false for text writing
 */
int pgmwrite(
   const char* filename,
   int w, int h, int pitch,
   unsigned char* data, 
   const char* comment_string = 0,
   bool binsave = true
);

/*!
 * \brief Write a PGM image from normalized floating point data.
 * 
 * The floats will be interpreted as black==0 and white==1,
 * values between are quantized to 8-bit grayscale, and
 * values outside that range will be clamped to black or white.
 * 
 * \param filename The file to write to.
 * \param w Image width
 * \param h Image height
 * \param data Row-major image data
 * \param comment_string Comments (NULL if none)
 * \param binsave 1 for binary writing, 0 for text writing
 */
int pgmwrite_float(
   const char* filename,
   int w, int h,
   float* data, 
   const char* comment_string,
   int binsave
);

/*!
 * \brief Write a PPM image.
 *
 * Writes binary format "P6".
 *
 * \param filename The file to write to.
 * \param w image width
 * \param h image height
 * \param pitch number of bytes in one row of data
 * \param data row-major image data
 * \param comment_string comments (NULL if none)
 */
int ppmwrite(
   const char* filename,
   int w, int h, int pitch,
   unsigned char const* data,
   const char* comment_string = 0
);

#endif /*PGVLPPM_H*/