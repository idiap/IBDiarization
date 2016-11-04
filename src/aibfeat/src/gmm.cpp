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


#include "gmm.h"

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
   float log_lkld= MINVALUEFORMINUSLOG ;
   for(int i=0; i<M; i++)
   {
      if(Weight[i])
      {
         local_lkld[i] = log(Weight[i]) + gArray[i].Log_Likelihood(feature);
         if(isnan(local_lkld[i]) || !finite(local_lkld[i]))
         {
            local_lkld[i] = MINVALUEFORMINUSLOG;
         }
         log_lkld = Log_Add(log_lkld, local_lkld[i]);
      }
      else{
         local_lkld[i] = MINVALUEFORMINUSLOG;
      }
   }
   return log_lkld;
}

void GMM::Initialize(float **data, int count)
{
   int gap = count/M;

   float *Mean = NULL;
   float *Var  = NULL;
   Mean = (float *)malloc(DIM*sizeof(float));
   Var  = (float *)malloc(DIM*sizeof(float));

   if((Mean == NULL) || (Var == NULL))
   {
      fprintf(stderr, "cannot allocate memory GMM::Initialize %lx %lx \n", 
            (unsigned long int)Mean, (unsigned long int)Var);
   }

   for(int i=0;i<M;i++)
   {
      Initialize1DFloat(Mean, DIM);
      Initialize1DFloat(Var, DIM);

      for(int k=0;k<gap;k++)
      {
         for(int j=0;j<DIM;j++)
         {
            Mean[j] += data[i*gap+k][j];
            Var[j] += pow(data[i*gap+k][j], 2.0);
         }  
      }

      for(int j=0;j<DIM;j++)
      {
         Mean[j] /= (float)gap;
         Var[j] /= (float)gap;
         Var[j] -= pow(Mean[j], 2.0);
         if(Var[j] <= MIN_POS_FLOAT)
         {
            Var[j] = MIN_POS_FLOAT;
         }
      }
      gArray[i].SetMean(Mean);
      gArray[i].SetVar(Var);
   }
   printf("Initializing the Gmm");
   free (Mean);
   free (Var);

}
float GMM::Train(float **data, int nSamples)
{
   float *Mean = NULL;
   float *Var  = NULL;
   Mean = (float *)malloc(DIM*sizeof(float));
   Var  = (float *)malloc(DIM*sizeof(float));

   if((Mean == NULL) || (Var == NULL))
   {
      fprintf(stderr, "cannot allocate memory GMM::Initialize %lx %lx \n", 
            (unsigned long int)Mean, (unsigned long int)Var);
   }

   if(M ==1)
   {
      Initialize1DFloat(Mean, DIM);
      Initialize1DFloat(Var, DIM);
      for(int i=0;i<nSamples;i++)
      {
         for(int j=0;j<DIM;j++)
         {
            Mean[j] += data[i][j];
            Var[j] += pow(data[i][j], 2.0);
         }
      }   
      for(int j=0;j<DIM;j++)
      {
         Mean[j] /= (float)nSamples;
         Var[j] /= (float)nSamples;
         Var[j] -= pow(Mean[j], 2.0);
         if(Var[j] <= MIN_POS_FLOAT)
         {
            Var[j] = MIN_POS_FLOAT;
         }
      }
      float det = 0.0;
      for(int i=0; i< DIM;i++)
      {
         det += log(Var[i]);
      }
      gArray[0].SetMean(Mean);
      gArray[0].SetVar (Var );
      free(Mean);
      free(Var);
      return -0.5*(float)nSamples*(det + DIM);
   }
   float **gamma_l_o = (float **)malloc(sizeof(float *)*M);
   float **gamma_l_o_2 = (float **)malloc(sizeof(float *)*M);
   float *gamma_l = (float *)malloc(sizeof(float)*M);       
   float total_log_lkld, lkld;
   float contri;

   *gamma_l_o = (float *)calloc(M*DIM, sizeof(float));
   *gamma_l_o_2 = (float *)calloc(M*DIM, sizeof(float));
   for(int j=1; j<M; j++)
   {
      gamma_l_o[j] = gamma_l_o[j-1] + DIM;
      gamma_l_o_2[j] = gamma_l_o_2[j-1] + DIM;
   }
   for(int loop=0; loop < N_LOOPS_ADAPT_GMM; loop++)
   {
      Initialize1DFloat(gamma_l, M);
      Initialize2DFloat(gamma_l_o, M, DIM);
      Initialize2DFloat(gamma_l_o_2, M, DIM);
      total_log_lkld = 0.0;
      for(int frame=0; frame<nSamples; frame++)
      {
         lkld = Log_Likelihood(data[frame]);
         total_log_lkld += lkld;
         for(int j=0; j<M; j++)
         {
            contri = local_lkld[j] - lkld;
            contri = (float)exp(contri);
            gamma_l[j] += contri;
            for(int fea=0; fea<DIM; fea++)
            {
               gamma_l_o[j][fea] += contri*data[frame][fea];
               gamma_l_o_2[j][fea] += contri*pow(data[frame][fea], 2.0);
            }
         }/*for(int j=0; j<M; j++)*/
      }/*for(int frame=0; frame<nSamples; frame++)*/
      for(int j=0; j<M; j++)
      {
         Weight[j] = gamma_l[j]/(float)nSamples;
         if(Weight[j] != 0.0)
         {
            for(int fea=0; fea <DIM; fea++)
            {
               Mean[fea] = gamma_l_o[j][fea]/gamma_l[j];
               Var[fea] = gamma_l_o_2[j][fea]/gamma_l[j];
               Var[fea] -= (float)pow(Mean[fea], 2.0);
               if(Var[fea] <= MIN_POS_FLOAT)
               {
                  Var[fea] = MIN_POS_FLOAT;
               }
            }
         }
         gArray[j].SetMean(Mean);
         gArray[j].SetVar(Var);

      }
   }/*for(int loop=0; loop < N_LOOPS_ADAPT; loop++)*/
   free(*gamma_l_o);free(*gamma_l_o_2);free(gamma_l_o);free(gamma_l_o_2);
   if(isnan(total_log_lkld) || !finite(total_log_lkld))
   {
      printf("The total log likelihood inside the GMM training invalid\n");
   }

   free (Mean);
   free (Var);
   return total_log_lkld;

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
   {
      gArray[i].SetMean(gmm.gArray[i].GetMean());
      gArray[i].SetVar (gmm.gArray[i].GetVar());
   }
}

GMM::GMM(const GMM &gmm)
{
   Initialize(gmm.M, gmm.DIM);
   memcpy(Weight, gmm.Weight, sizeof(float)*M);
   lkld = gmm.lkld;
   for(int i=0; i<M; i++)
   {
      gArray[i].SetMean(gmm.gArray[i].GetMean());
      gArray[i].SetVar (gmm.gArray[i].GetVar());
   }
}

float GMM::Log_Add(float log_a, float log_b)
{
   float result;
   if(log_a < log_b)
   {
      float tmp = log_a;
      log_a = log_b;
      log_b = tmp;
   }
   if((log_b - log_a) <= MINVALUEFORMINUSLOG)
   {
      return log_a;
   }
   else
   {
      result = log_a + (float)log(1.0 + exp(log_b - log_a));
   }
   return result;
}

void Gaussian::Initialize(int dim)
{
   DIM = dim;
   if(Mean != NULL)
      delete Mean;
   if(Var != NULL)
      delete []Var;
   if(logVar != NULL)
      delete logVar;

   Mean = new float[DIM];
   Var = new float[DIM];
   logVar = new float[DIM];

   //printf("The Gaussian component is initialized\n");

}
void Gaussian::SetMean(const float* mean)
{
   memcpy(Mean,mean, sizeof(float)*DIM);
}

void Gaussian::SetVar(const float* var)
{
   int i;
   memcpy(Var,var, sizeof(float)*DIM);
   for(i=0; i<DIM; ++i)
      logVar[i] = log(Var[i]);
}
const float* Gaussian::GetMean()
{
   return Mean;
}
const float* Gaussian::GetVar()
{
   return Var;
}
Gaussian::Gaussian()
{
   Mean = NULL;
   Var  = NULL;
   logVar = NULL;
   DIM = 0;
}

Gaussian::~Gaussian()
{
   if(Mean != NULL)
      delete []Mean;
   if(Var != NULL)
      delete []Var;
}
float Gaussian::Log_Likelihood(float *feature)
{
   float x,y=0;
   for(int i=0; i<DIM; i++)
   {
      x = feature[i]-Mean[i];
      y += x*x/Var[i]+logVar[i];
   }
   return((float)-0.5*y);
}

