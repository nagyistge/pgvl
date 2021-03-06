/* ppm.cpp is part of pgvl and is 
 * Copyright 2015 Philip G. Lee <rocketman768@gmail.com>
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

unsigned char* pgmread(const char* filename, int* w, int* h)
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
    
    if ((data = new unsigned char[numpix]()) == NULL)
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

float* pgmread_float(const char* filename, int* w, int* h )
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

unsigned char* ppmread(const char* filename, int* w, int* h, int* maxval)
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
    
    if ((data = new unsigned char[numpix*3]()) == NULL)
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

float* ppmread_float(const char* filename, int* w, int* h )
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

int pgmwrite(
   const char* filename,
   int w, int h, int pitch,
   unsigned char* data, 
   const char* comment_string,
   bool binsave
)
{
    FILE* file;
    int maxval;
    int nread;
    int i,j,k;
    
    if ((file = fopen(filename, "w")) == NULL)
    {
       printf("ERROR: file open failed\n");
       return(-1);
    }

    if (binsave)
      fprintf(file,"P5\n");
    else
      fprintf(file,"P2\n");

    if (comment_string)
      fprintf(file,"# %s\n", comment_string);

    fprintf(file,"%d %d\n", w, h);
    
    maxval = 255;
    fprintf(file, "%d\n", maxval);
    
    if (binsave)
    {
      while(h--)
      {
        nread = fwrite(data, sizeof(unsigned char), w, file);
        if( nread != w )
        {
          fprintf(stderr, "Error: wrote %d/%d pixels.", nread, w);
          exit(1);
        }

        data += pitch;
      }
    }
    else
    {
      printf("Writing to %s as ascii.\n", filename);

      for(i=0; i<h; i++)
      {
        k = 0;
        for(j=0; j<w; j++)
          fprintf(file,"%d ", (int)data[k++]);
        data += pitch;
      }
    }     
   
    fclose(file);
    return 0;
}

int pgmwrite_float(
   const char* filename,
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
   ret = pgmwrite(filename, w, h, w, cdata, comment_string, binsave);
   
   delete[] cdata;
   return ret;
}

int ppmwrite(
   const char* filename,
   int w, int h, int pitch,
   unsigned char const* data,
   const char* comment_string
)
{
    FILE* file;
    int maxval;
    int nread;
    int rowpix = 3*w;

    if ((file = fopen(filename, "w")) == NULL)
    {
       printf("ERROR: file open failed\n");
       return(-1);
    }

    fprintf(file,"P6\n");
    if (comment_string)
      fprintf(file,"# %s\n", comment_string);
    fprintf(file,"%d %d\n", w, h);

    maxval = 255;
    fprintf(file, "%d\n", maxval);

    while(h--){
      nread = fwrite(data, sizeof(unsigned char), rowpix, file);
      if( nread != rowpix )
      {
        fprintf(stderr, "Error: wrote %d/%d pixels.", nread, rowpix);
        exit(1);
      }
      data += pitch;
    }

    fclose(file);
    return 0;
}
