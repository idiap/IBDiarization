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
#include "rgmm.h"
#include <vector>

#define LINE_LENGTH 1024

using namespace std;

namespace gaussrealignment
{


GMM::GMM()
{
}

void GMM::Initialize(int num, int dim)
{
   M = num; 
   DIM = dim;
   Weight = new float[M];
   gArray = new Gaussian[M];
   local_lkld = new float[M];
   for(int i=0; i<M; i++)
   {
      Weight[i] = 1.0/(float)M;
      gArray[i].Initialize(dim);
   }
}

float GMM::Log_Likelihood(float *feature)
{
   return gArray[0].Log_Likelihood(feature); //KL HMM has only one state distribution
}

void GMM::Initialize(float **data, int count)
{
   //We make use of the fact that number of mixtures is 1
   //

   gArray[0].Initialize(DIM);

   float *mean = new float[DIM];

   for(int i=0;i<count;i++)
   {
      for(int j=0;j<DIM;j++)
         mean[j] += data[i][j];
   }

   for(int j=0;j<DIM;j++)
      mean[j] /= count;

   gArray[0].SetMean(mean);

}



float GMM::Train(float **data, int nSamples)
{
   // The new training has to train only one mean vector 
   // per GMM ... (GMM has one gaussian and it contains the 
   // mean posterior distribution 
   int i;
   float *Mean = NULL; 
   Mean = (float *)malloc(DIM*sizeof(float));
   
   

   if(Mean == NULL)
   {
      fprintf(stderr, "cannot allocate memory GMM::Train for Mean\n");
      exit(0);
   }

   for (i = 0; i <DIM; ++i)
      Mean[i] = 0;
   // Initialize the models on the data
   //Initialize(data,nSamples); seemed like a good idea, but scores are worse.
   for(int frame=0; frame<nSamples; frame++) {
      float *curr_data = data[frame];
      for (i = 0; i<DIM; ++i)
      {
         Mean[i] += curr_data[i];
         //Var[i] += (curr_data[i]* curr_data[i]);
      }
   }/*for(int frame=0; frame<nSamples; frame++)*/

   float sum = 0;

   for (i = 0; i<DIM; ++i)
   { 
      if(Mean[i] < 0)
         Mean[i] = 0;

      sum += Mean[i];
   }

   for (i = 0; i<DIM; ++i) {
      Mean[i] /= sum;
   }
   gArray[0].SetMean(Mean);
   //gArray[0].SetVar(Var);
   // Evaluate the newly trained model on the data
   //

   free(Mean);
   return 0.0;
}

GMM::~GMM()
{
   delete []gArray;
   delete []Weight;
   delete []local_lkld;
}

void GMM::Copy(GMM gmm)
{
   Initialize(gmm.M, gmm.DIM);
   memcpy(Weight, gmm.Weight, sizeof(float)*M);
   lkld = gmm.lkld;
   for(int i=0; i<M; i++)
      gArray[i].SetMean(gmm.gArray[i].GetMean());
}

GMM::GMM(const GMM &gmm)
{
   Initialize(gmm.M, gmm.DIM);
   memcpy(Weight, gmm.Weight, sizeof(float)*M);
   lkld = gmm.lkld;
   for(int i=0; i<M; i++)
      gArray[i].SetMean(gmm.gArray[i].GetMean());
}

float GMM::Log_Add(float log_a, float log_b)
{
   float result;
   if(log_a < log_b) {
      float tmp = log_a;
      log_a = log_b;
      log_b = tmp;
   }
   if((log_b - log_a) <= MINVALUEFORMINUSLOG) {
      return log_a;
   }
   else {
      result = log_a + (float)log(1.0 + exp(log_b - log_a));
   }
   return result;
}

void Gaussian::SetMean(const float* mean)
{
   int i;
   memcpy(Mean, mean, sizeof(float)*DIM);
   for(i=0; i<DIM; ++i) {
      logMean[i] = log(Mean[i]);
      if(isnan(logMean[i])) {
          cout << "Warning: Mean is nan!" << endl;
      }
      else if(isinf(logMean[i])) {
          cout << "Warning: Mean is inf!" << endl;
      }
   }
}

const float* Gaussian::GetMean()
{
   return Mean;
}
const float* Gaussian::GetLogMean()
{
   return logMean;
}

void Gaussian::Initialize(int dim)
{
   DIM = dim;
   Mean = new float[DIM];
   logMean = new float[DIM];
}

Gaussian::Gaussian()
{
   Mean = NULL;
   logMean = NULL;
}

Gaussian::~Gaussian()
{
   if(Mean != NULL) delete []Mean;
   if(logMean != NULL) delete []logMean;
}
float Gaussian::Log_Likelihood(float *feature)
{
   float y = 0.0;
   for(int i=0; i<DIM; i++)
   {
      if (feature[i] > 0 && Mean[i] > 0)
      {
         y += (feature[i] * (-logMean[i]));
      }
   }
   if(isnan(y) || isinf(y)) {
       cout << "Warning: error in likelihood value = " << y << endl;
       y = 10E8;        //set a large value
   }
   return((float)-y); //Multiplying with \beta
}

int realign_read_scp_file(const char *in_filename, int **out_seginfo, int *num_segs)
{
   char buffer[LINE_LENGTH];
   char tmp[LINE_LENGTH];
   int start,end;
   vector< vector<int> > pos;
   vector<int> ivel;

   FILE* scpfp = NULL;
   scpfp = fopen(in_filename,"r");
   if(scpfp == NULL)
   {
      fprintf(stderr,"Cannot Open the scp file<%s>\n",in_filename);
      return -1;
   }
   fgets(buffer,LINE_LENGTH,scpfp);
   while(!feof(scpfp))
   {
      ivel.clear();
      sscanf(buffer,"%[^[][%d,%d",tmp,&start,&end);
      fgets(buffer,LINE_LENGTH,scpfp);
      ivel.push_back(start);
      ivel.push_back(end);
      pos.push_back(ivel);
   }

   *out_seginfo = (int*) malloc(2*pos.size()*sizeof(int));
   for (unsigned int i =0; i <pos.size(); ++i)
   {
      (*out_seginfo)[2*i] = pos[i][0];
      (*out_seginfo)[2*i+1] = pos[i][1];
   }

   *num_segs = pos.size();
   fclose(scpfp);
   return 0;
}
}
