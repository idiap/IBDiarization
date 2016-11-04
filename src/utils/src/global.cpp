// Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
// Written by Fabio Valente <fabio.valente@idiap.ch>
// Written by Deepu Vijayasenan <dvijayasenan@lsv.uni-saarland.de>
// Written by David Imseng <david.imseng@idiap.ch>
// Written by Srikanth Madikeri <srikanth.madikeri@idiap.ch>
//
// This file is part of the IB Speaker Diarization Toolkit.
//
// IB diarization toolkit is free software: you can redistribute it
// and/or modify it under the terms of the GNU General Public License
// version 3 as published by the Free Software Foundation.
//
// The IB Speaker Diarization Toolkit is distributed in the hope that it
// will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the IB Speaker Diarization Toolkit. If not, see
// <http://www.gnu.org/licenses/>.


#include "global.h"

void initializeF(float **mem, int size)
{
   *mem = NULL;
   *mem = (float *) malloc (size * sizeof(float));
   if (*mem == NULL )
      fprintf(stderr, "ERROR in float memory init");
}

void initializeI(int **mem, int size)
{
   *mem = NULL;
   *mem = (int *) malloc (size * sizeof(int));
   if (*mem == NULL )
      fprintf(stderr, "ERROR in int memory init");
}

void releaseF(float **mem)
{
   if (*mem != NULL)
      free (*mem);
   *mem = NULL;
}

void releaseI(int **mem)
{
   if (*mem != NULL)
      free (*mem);
   *mem = NULL;
}

void Initialize1DFloat(float *array, int dim)
{
   for(int i=0; i<dim; i++)
   {
      array[i] = 0.0;
   }
}

void Initialize1DInt(int *array, int dim)
{
   for(int i=0; i<dim; i++)
   {
      array[i] = 0;
   }
}

void Initialize2DFloat(float **array, int dim1, int dim2)
{
   for(int i=0; i<dim1; i++)
   {
      for(int j=0; j<dim2; j++)
      {
         array[i][j] = 0.0;
      }
   }
}


void subtract(float *array1, float *array2, int dim)
{
   for(int i=0;i<dim;i++)
   {
      array1[i] -= array2[i];
   }
}

int SwapByteOrderOfInt(void *Source)
{
   // These definitions assume unsigned char is 8 bits wide
   // and that each entry in an array of them is at the byte
   // address immediately following the previous entry
   unsigned char *From = (unsigned char *)Source;
   int RetVal;
   unsigned i;
   unsigned char *To = (unsigned char *)&RetVal;

   // For each byte of the source
   for (i = 0; i < sizeof(int); i++)
   {
      // Copy each byte to the target in reverse order
      To[sizeof(int) - i - 1] = From[i];
   }
   return RetVal;
}


__uint32_t SwapByteOrderOfLong(void *Source)
{
   // These definitions assume unsigned char is 8 bits wide
   // and that each entry in an array of them is at the byte
   // address immediately following the previous entry
   unsigned char *From = (unsigned char *)Source;
   __uint32_t RetVal;
   unsigned i;
   unsigned char *To = (unsigned char *)&RetVal;

   // For each byte of the source
   for (i = 0; i < SIZE_LONG ; i++)
   {
      // Copy each byte to the target in reverse order
      To[SIZE_LONG  - i - 1] = From[i];
   }
   return RetVal;
}

float SwapByteOrderOfFloat(void *Source)
{
   // These definitions assume unsigned char is 8 bits wide
   // and that each entry in an array of them is at the byte
   // address immediately following the previous entry
   unsigned char *From = (unsigned char *)Source;
   float RetVal;
   unsigned i;
   unsigned char *To = (unsigned char *)&RetVal;

   // For each byte of the source
   for (i = 0; i < sizeof(float); i++)
   {
      // Copy each byte to the target in reverse order
      To[sizeof(float) - i - 1] = From[i];
   }
   return RetVal;
}


__uint16_t SwapByteOrderOfShort(void *Source)
{
   // These definitions assume unsigned char is 8 bits wide
   // and that each entry in an array of them is at the byte
   // address immediately following the previous entry
   unsigned char *From = (unsigned char *)Source;
   __uint16_t RetVal;
   unsigned i;
   unsigned char *To = (unsigned char *)&RetVal;

   // For each byte of the source
   for (i = 0; i < SIZE_SHORT; i++)
   {
      // Copy each byte to the target in reverse order
      To[SIZE_SHORT - i - 1] = From[i];
   }
   return RetVal;
}

