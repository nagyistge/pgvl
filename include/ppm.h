/* ppm.h is part of pgvl and is 
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */

#include <string.h>
#include <stdio.h>

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
unsigned char* pgmread(char* filename, int* w, int* h)
{
    FILE* file;
    char line[256];
    int maxval;
    int binary;
    int nread;
    int numpix;
    int i,j,k,int_tmp;

    unsigned char* data;
    
    if ((file = fopen(filename, "r")) == NULL)
    {
       printf("ERROR: file open failed\n");
       *h = *w = 0;
       return(NULL);
    }
    if( fgets(line, 256, file) != line )
    {
       fprintf(stderr, "ERROR: something wrong with the file.\n");
       return NULL;
    }
    if (strncmp(line,"P5", 2))
    {
       if (strncmp(line,"P2", 2))
       {
          printf("pgm read: not a pgm file\n");
          *h = *w = 0;
          return(NULL);
       }
       else 
          binary = 0;
    }
    else 
       binary = 1;

    if( fgets(line, 256, file) != line )
    {
       fprintf(stderr, "ERROR: something wrong with the file.\n");
       return NULL;
    }
    while (line[0] == '#')
    {
       if( fgets(line, 256, file) != line )
       {
          fprintf(stderr, "ERROR: something wrong with the file.\n");
          return NULL;
       }
    }

    sscanf(line,"%d %d", w, h);
    if( fgets(line, 256, file) != line )
    {
       fprintf(stderr, "ERROR: something wrong with the file.\n");
       return NULL;
    }
    sscanf(line, "%d", &maxval);
    
    numpix = (*w)*(*h);
    
    if ((data = (unsigned char*)calloc(numpix, sizeof(unsigned char*))) == NULL)
    {
      printf("Memory allocation error. Exit program");
      exit(1);
    }

    if (binary)
    {
       nread = fread( (void*)data, sizeof(unsigned char), numpix, file);
       if( nread != numpix )
       {
         fprintf(stderr, "Error: read %d/%d pixels.", nread, numpix);
         exit(1);
       }
    }
    else
    {
       k = 0;
       for (i = 0; i < (*h); i++)
       {
          for (j = 0; j < (*w); j++)
          {
             if( fscanf(file, "%d", &int_tmp) != 1 )
             {
                fprintf(stderr, "ERROR: Something wrong with the file.\n");
                exit(1);
             }
             data[k++] = (unsigned char)int_tmp;
          }  
       }
    }
    
    fclose(file);
    return data;
}

/*!
 * \brief Read a normalized grayscale floating-point image.
 * 
 * \sa pgmread()
 */
float* pgmread_float(char* filename, int* w, int* h )
{
   int i, numpix;
   unsigned char* cdata;
   float* fdata;
   
   cdata = pgmread(filename, w, h);
   numpix = (*w)*(*h);
   fdata = new float[numpix];
   
   for( i = 0; i < numpix; ++i )
      fdata[i] = static_cast<float>(cdata[i])/255.0f;
   
   delete[] cdata;
   return fdata;
}

unsigned char* rgbToArgb(unsigned char* rgb, int w, int h)
{
   int j,k;
   unsigned char* ret = new unsigned char[4*w*h];
   
   for( k=0,j=0; k < 4*w*h; ++k )
   {
      if( k%4 == 0 )
         ret[k] = 0xff;
      else
         ret[k] = rgb[j++];
   }
   
   return ret;
}

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
unsigned char* ppmread(char* filename, int* w, int* h, int* maxval)
{
    FILE* file;
    char line[256];
    int binary;
    int nread;
    int numpix;
    int i,j,k,int_tmp;

    unsigned char* data;
    
    if ((file = fopen(filename, "r")) == NULL)
    {
       printf("ERROR: file open failed\n");
       *h = *w = 0;
       return(NULL);
    }
    if( fgets(line, 256, file) != line )
    {
       fprintf(stderr, "ERROR: something wrong with the file.\n");
       return NULL;
    }
    if (strncmp(line,"P6", 2))
    {
       if (strncmp(line,"P3", 2))
       {
          printf("ppm read: not a ppm file\n");
          *h = *w = 0;
          return(NULL);
       }
       else 
          binary = 0;
    }
    else 
       binary = 1;

    if( fgets(line, 256, file) != line )
    {
       fprintf(stderr, "ERROR: something wrong with the file.\n");
       return NULL;
    }
    while (line[0] == '#')
    {
       if( fgets(line, 256, file) != line )
       {
          fprintf(stderr, "ERROR: something wrong with the file.\n");
          return NULL;
       }
    }

    sscanf(line,"%d %d", w, h);
    if( fgets(line, 256, file) != line )
    {
       fprintf(stderr, "ERROR: something wrong with the file.\n");
       return NULL;
    }
    sscanf(line, "%d", maxval);
    
    if( *maxval < 0 || *maxval > 255 )
    {
       fprintf(stderr, "Error: maximum value %d is bad.\n", *maxval);
       exit(1);
    }
    
    numpix = (*w)*(*h);
    
    if ((data = (unsigned char*)calloc(numpix*3, sizeof(unsigned char*))) == NULL)
    {
      printf("Memory allocation error. Exit program");
      exit(1);
    }

    if (binary)
    {
       nread = fread( (void*)data, sizeof(unsigned char), numpix*3, file);
       if( nread != numpix*3 )
       {
          fprintf(stderr, "Error: read %d/%d pixels.", nread/3, numpix);
          exit(1);
       }
    }
    else
    {
       k = 0;
       for (i = 0; i < (*h); i++)
       {
          for (j = 0; j < (*w)*3; j++)
          {
             if( fscanf(file, "%d", &int_tmp) != 1 )
             {
                fprintf(stderr, "ERROR: Something wrong with the file.\n");
                exit(1);
             }
             data[k++] = (unsigned char)int_tmp;
          }  
       }
    }
    
    fclose(file);
    return data;
}

/*!
 * \brief Read a normalized floating-point image.
 * 
 * \sa ppmread()
 */
float* ppmread_float(char* filename, int* w, int* h )
{
   int maxval;
   int i, numpix;
   unsigned char* cdata;
   float* fdata;
   
   cdata = ppmread(filename, w, h, &maxval);
   numpix = (*w)*(*h);
   fdata = new float[3*numpix];
   
   for( i = 0; i < numpix*3; ++i )
      fdata[i] = static_cast<float>(cdata[i])/maxval;
   
   delete[] cdata;
   return fdata;
}

/*!
 * \brief Write a PGM image.
 * 
 * \param filename The file to write to.
 * \param w Image width
 * \param h Image height
 * \param data Row-major image data
 * \param comment_string Comments (NULL if none)
 * \param binsave 1 for binary writing, 0 for text writing
 */
int pgmwrite(
   char* filename,
   int w, int h,
   unsigned char* data, 
   const char* comment_string,
   int binsave
)
{
    FILE* file;
    int maxval;
    int nread;
    int i,j,k;
    int numpix = w*h;
    
    if ((file = fopen(filename, "w")) == NULL)
    {
       printf("ERROR: file open failed\n");
       return(-1);
    }

    if (binsave == 1)
      fprintf(file,"P5\n");
    else
      fprintf(file,"P2\n");

    if (comment_string != NULL)
      fprintf(file,"# %s \n", comment_string);

    fprintf(file,"%d %d \n", w, h);
    
    // I think this is wrong.
    /*
    maxval = 0;
    k = 0;
    for (i = 0; i < h; i++)
    {
       for (j=0; j < w; j++)
       {
         if ((int)data[k] > maxval)
            maxval = (int)data[k];
         ++k;
       }
    }
   */
    
    maxval = 255;
    fprintf(file, "%d \n", maxval);
    
    if (binsave == 1)
    {
      nread = fwrite(data, sizeof(unsigned char), numpix, file);
      if( nread != numpix )
      {
         fprintf(stderr, "Error: wrote %d/%d pixels.", nread, numpix);
         exit(1);
      }
    }
    else
    {
      printf("Writing to %s as ascii.\n", filename);

      k = 0;
      for(i=0; i<h; i++)
        for(j=0; j<w; j++)
          fprintf(file,"%d ", (int)data[k++]);
    }     
   
    fclose(file);
    return 0;
}

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
   char* filename,
   int w, int h,
   float* data, 
   const char* comment_string,
   int binsave
)
{
   int i, numpix;
   int ret;
   unsigned char* cdata;
   
   numpix = w*h;
   cdata = new unsigned char[numpix];
   
   for( i = 0; i < numpix; ++i )
      cdata[i] = data[i] < 0.f ? 0 : (data[i]>1.f?255:(static_cast<unsigned char>(255.f*data[i])));
   ret = pgmwrite(filename, w, h, cdata, comment_string, binsave);
   
   delete[] cdata;
   return ret;
}