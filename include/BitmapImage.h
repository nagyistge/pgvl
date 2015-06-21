/*
 * BitmapImage.h is part of pgvl and is
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */

#ifndef BITMAPIMAGE_H
#define BITMAPIMAGE_H

#include <string>
#include <inttypes.h>

/*!
 * \brief An image class for bitmaps
 *
 * For efficiency, each row is aligned to \c CACHE_LINE_SIZE byte boundaries.
 */
class BitmapImage {
public:
   //! \brief Default constructor
   BitmapImage();
   ~BitmapImage();

   /*!
    * \brief Constructor from image filename
    *
    * \param filename a .ppm or a .pgm image
    */
   BitmapImage(std::string const& filename);

   //! \brief Number of rows in the image
   int rows() const { return _rows; }
   //! \brief Number of columns in the image
   int cols() const { return _cols; }
   //! \brief Number of channels in the image
   int channels() const { return _channels; }
   //! \brief Number of bytes in a row of pixels
   int rowWidth() const { return _rowWidth; }

   /*!
    * \brief Save the file
    *
    * If the number of channels is 1, it appends a .pgm extension and writes.
    * If the number of channels is 3, it appends a .ppm extension and writes.
    *
    * \param basename the base filename without extension
    */
   void save(std::string const& basename) const;

   //! \brief Pointer to the ith row of pixel data
   uint8_t* operator[](size_t i) { return _data + i*_rowWidth; }
   //! \brief Pointer to the ith row of pixel data (const version)
   uint8_t const* operator[](size_t i) const { return _data + i*_rowWidth; }
private:

   uint8_t* _data;
   uint8_t* _unalignedData;
   int _rows;
   int _cols;
   int _channels;
   int _rowWidth;

   enum FileType { FILETYPE_NONE, FILETYPE_PGM, FILETYPE_PPM };
   static FileType fileType(std::string const& filename);
};

#endif /*BITMAPIMAGE_H*/
