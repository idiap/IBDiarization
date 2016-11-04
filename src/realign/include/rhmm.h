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
#include "featfileread.h"
#include "rgmm.h"

using namespace gaussrealignment;

namespace hmmrealignment
{

class chunk
{
   public:
      int from;
      int to;
      int cluster;
      chunk *next;
      chunk *prev;
      ~chunk()
      {
         if(next!= NULL)
            delete next;
      }
};

class HMM
{
   public:
      int nClass;
      int DIM;
      int MIN_DURATION;
      float **Clusters_Means;
      float *lkld;
      GMM *gmm;
      int **NP;
      int *U; /* a pointer to indicate if a class is usable or not*/
      HMM(int K, int M, int DIM, int MD);
      ~HMM();
      HMM();
      float segment(FeatFileReader& f, chunk **pFirst);
      void Calculate_Lkld(float *data);
      void CMS(FeatFileReader& f, chunk *);
      void Train_HMM(FeatFileReader& f, chunk *first);

      void Get_Data(FeatFileReader& f, float ***pData, chunk *first, int cd1, int cd2, int *count);
      int Find_N_Best(int *mcd1, int *mcd2, int *N, chunk *first, char *fname);
      void Write(char *fname, chunk *first);
      void WriteLkld(char *fname, char *name_dat);
      void WriteParameters(char *);
      int if_possible(int cd1, int cd2);
      void WriteRttm(const char *fname, chunk *first, FeatFileReader& fileReader, int tot_num_frames, const char *mtgid);
};


}
