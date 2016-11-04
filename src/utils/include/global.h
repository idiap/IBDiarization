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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

//#define MIN_DURATION 200 // number of frames
#define LOG_0_5 -0.693147181
#define LOG_0_9 -0.105360516//-0.693147181//-0.510825624//-0.223143551//-0.356674944//-0.105360516
#define LOG_0_1 -2.302585093//-0.693147181//-0.916290732//-1.609437912//-1.203972804//-2.302585093
#define NLOOPS_ADAPT 1
#define N_LOOPS_ADAPT_GMM 5
#define MINVALUEFORMINUSLOG -103.0 // -18.42
#define MIN_POS_FLOAT 1.4e-45
#define N_LOOPS_K_MEANS 15
#define LAST 0
#define CURR 1
#define P 2	//variable to decide the initialization
#define MIN_FLOAT -3.40282346638528860e+38
#define FULL_HEADER 1
#define ONLY_SEGMENTATION 0
#define K_MEANS 0
#define SIZE_LONG 4
#define SIZE_SHORT 2
#define VAR_TH 1e-06

void initializeF(float **mem,int size);
void initializeI(int **mem,int size);

void releaseF(float **mem);
void releaseI(int **mem);

void Initialize1DFloat(float *array, int dim);
void Initialize1DInt(int *array, int dim);
void Initialize2DFloat(float **array, int dim1, int dim2);
void subtract(float *, float *, int);

//Internal functions for Byte Swapping
__uint16_t SwapByteOrderOfShort(void *Source);
int SwapByteOrderOfInt(void *Source);
__uint32_t SwapByteOrderOfLong(void *Source);
float SwapByteOrderOfFloat(void *Source);

