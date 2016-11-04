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

class FeatFileReader {

   FILE *m_featfp; 
   int *m_indexMap;
   __uint16_t m_byte_per_frame;
   int m_featdim;
   int m_headersize;
   int m_nSamples;
   int m_sampleRate;

   public:
   FeatFileReader();
   FeatFileReader(const char* featfilename, const int* seg_info, int num_segs, int& featdim, int& num_fr);
   ~FeatFileReader();
   void readfeat(int fr_index,float *featvec);
   int get_nSamples();
   int getSampleRate();
   int findFrameIndex(int ); //This function accepts the speech only frame index as input and return the actual frame index
};
