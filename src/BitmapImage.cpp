#include <BitmapImage.h>
#include <ppm.h>
#include "config.h"

#include <regex>
#include <iostream>
#include <stdlib.h>

BitmapImage::BitmapImage() :
   _data(0),
   _unalignedData(0),
   _rows(0),
   _cols(0),
   _channels(0),
   _rowWidth(0)
{
}

BitmapImage::~BitmapImage() {
   delete[] _unalignedData;
}

BitmapImage::BitmapImage(std::string const& filename) :
   _data(0),
   _unalignedData(0),
   _rows(0),
   _cols(0),
   _channels(0),
   _rowWidth(0)
{
   int maxval;
   uint8_t* rawData = 0;

   switch( fileType(filename) ) {
   case FILETYPE_PGM:
      rawData = pgmread(filename.c_str(), &_rows, &_cols);
      _channels = 1;
      break;
   case FILETYPE_PPM:
      rawData = ppmread(filename.c_str(), &_rows, &_cols, &maxval);
      _channels = 3;
      break;
   default:
      LOGE("Bad image format");
      return;
   }

   // Make the row width a multiple of CACHE_LINE_SIZE
   _rowWidth = _channels*_cols;
   if( _rowWidth % CACHE_LINE_SIZE )
      _rowWidth += CACHE_LINE_SIZE - (_rowWidth % CACHE_LINE_SIZE);

   // Make each row line up with CACHE_LINE_SIZE
   _unalignedData = new uint8_t[_rowWidth*_rows + CACHE_LINE_SIZE];
   _data = _unalignedData;
   if( reinterpret_cast<size_t>(_data) % CACHE_LINE_SIZE )
      _data += CACHE_LINE_SIZE - (reinterpret_cast<size_t>(_data) % CACHE_LINE_SIZE);

   // Copy the data to aligned memory
   for( int i = 0; i < _rows; ++i ) {
      memcpy(_data + i*_rowWidth, rawData + i*_channels*_cols, _channels*_cols);
   }

   delete[] rawData;
}

void BitmapImage::save(std::string const& basename) const {
   uint8_t* rawData = new uint8_t[_rows*_cols*_channels];

   for( int i = 0; i < _rows; ++i ) {
      memcpy(rawData + i*_channels*_cols, _data + i*_rowWidth, _channels*_cols);
   }

   std::string filename(basename);

   switch( _channels ) {
   case 1:
      filename += ".pgm";
      pgmwrite(filename.c_str(), _cols, _rows, rawData, "", 1);
      break;
   case 3:
      // TODO: implement ppmwrite()
      //filename += ".ppm";
      //ppmwrite(filename.c_str(), _cols, _rows, rawData, "", 1);
      LOGE("PPM writing not yet supported");
      break;
   default:
      LOGE("Bad image format");
      break;
   }

   delete[] rawData;
}

BitmapImage::FileType BitmapImage::fileType(std::string const& filename) {
   std::regex ppm(".*[.]ppm$");
   std::regex pgm(".*[.]pgm$");

   if( std::regex_match(filename, pgm) )
      return FILETYPE_PGM;
   else if( std::regex_match(filename, ppm) )
      return FILETYPE_PPM;
   else
      return FILETYPE_NONE;
}