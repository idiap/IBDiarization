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
#include "rhmm.h"

using namespace gaussrealignment;

namespace hmmrealignment
{

HMM::HMM(int number_of_classes, int M, int dim, int md)
{
   nClass = number_of_classes;
   DIM = dim;
   MIN_DURATION = md;
   U = new int[nClass];
   gmm = new GMM[nClass];
   lkld = new float[nClass];
   Clusters_Means= (float **)malloc(sizeof(float *)*nClass);
   *Clusters_Means= (float *)malloc(sizeof(float)*nClass*DIM);
   NP = (int **)malloc(sizeof(int *)*nClass);
   *NP = (int *)malloc(sizeof(int)*nClass*nClass);
   for(int i=1;i<nClass;i++)
   {
      NP[i] = NP[i-1] + nClass;
      Clusters_Means[i] = Clusters_Means[i-1] + DIM;
   }
   Initialize2DFloat(Clusters_Means, nClass, DIM);
   for(int i=0; i<nClass; i++)
   {
      NP[i][0] = 0;//the first element is going to indicate the number
      gmm[i].Initialize(M, DIM);
      U[i] = 1;
   }
}

HMM::~HMM()
{
   delete []U;
   delete []gmm;
   delete []lkld;
   free(*NP);free(NP);
   free(*Clusters_Means);free(Clusters_Means);
}

HMM::HMM()
{
}


float HMM::segment(FeatFileReader &featFile, chunk **sent)
{
   //FILE *fin = fopen(fname, "rb");   
   int nSamples, MaxClass;
   float MaxLkld;
   int cl; //index variables cl = class index variable
   //st state index variable
   nSamples = featFile.get_nSamples();
   float *dC[2];
   int *tC[2];
   int totStates = (MIN_DURATION)*nClass;
   int *F = new int[nSamples+1];
   int *R = new int[nSamples+1];
   dC[LAST] = new float[totStates];
   dC[CURR] = new float[totStates]; 
   tC[LAST] = new int[totStates]; 
   tC[CURR] = new int[totStates]; 
   R[0] = -1; // cluster ids of frames
   F[0] = 1; //Frame indices for back tracking 
   int n_usable = 0;
   for(int cl=0; cl<nClass; cl++)
   {
      if(U[cl])
         n_usable++;
   }


   printf("For the Segmentation, nSamples = %d, nClass = %d\n", nSamples, n_usable);
   fflush(stdout);
   float TRANS_PROB = log((float)0.1/(float)n_usable);
   float *data = new float[DIM];
   int st=0; 
   featFile.readfeat(0,data);
   Calculate_Lkld(data);


   for(st=0; st<totStates; st++)
   {
      dC[LAST][st] = -MINVALUEFORMINUSLOG*totStates; //Just a large negative value
      dC[CURR][st] = -MINVALUEFORMINUSLOG*totStates; //Just a large negative value
      tC[LAST][st] = 0;
      tC[CURR][st] = 0;
   }

   for(st=0; st<nClass; st++)
   {
      dC[LAST][(st*(MIN_DURATION))] = (double)(lkld[st]);
      dC[CURR][(st*(MIN_DURATION))] = (double)(lkld[st]);
   }


   for(int frame = 1; frame <nSamples; frame++)
   {
      featFile.readfeat(frame,data);
      Calculate_Lkld(data); // now lkld[st] gives the values of each likelihood of this frame


      //State 0 of each class  is calculated from the 
      //max over all those previous states 

      // Step -  I calculate the max over the states... 
      double max_lkld_value  = 0;
      if (nClass>0)
      {
         max_lkld_value = dC[LAST][MIN_DURATION -1 ];
         for (cl = 1; cl<nClass; cl++)
         {
            //Comparing the last state likelihood values
            if (U[cl] && (max_lkld_value < dC[LAST][cl*MIN_DURATION + MIN_DURATION -1] ))
            {
               max_lkld_value = dC[LAST][cl*MIN_DURATION + MIN_DURATION -1 ];
            }
         }

      }
      // Step - II copy to all classes' first state with transition prob and likelihood attached
      for (cl=0; cl <nClass; ++cl)
      {
         st = cl*MIN_DURATION;
         dC[CURR][st] = max_lkld_value + TRANS_PROB + lkld[cl];
         tC[CURR][st] = frame; //Since we are transitioning this frame is the starting of a segment
      }



      // 2 - N-1 states of each class  (min duration )
      // likelihood = prev_lkld 
      for (cl=0; cl <nClass; ++cl)
      {
         for (st=1; st<MIN_DURATION-1; ++st)
         {
            int currStateIndx = cl*MIN_DURATION + st; 
            dC[CURR][currStateIndx] =  dC[LAST][currStateIndx-1]+ lkld[cl];
            tC[CURR][currStateIndx] = tC[LAST][currStateIndx-1]; //not transitioning, starting fr is same as before
         }
      }

      //Now the last state 
      for (cl=0; cl <nClass; ++cl)
      {
         int currStateIndx = cl*MIN_DURATION + MIN_DURATION -1; 
         if(dC[LAST][currStateIndx-1] > dC[LAST][currStateIndx] + LOG_0_9)
         {
            dC[CURR][currStateIndx] =  dC[LAST][currStateIndx-1]+ lkld[cl];
            tC[CURR][currStateIndx] = tC[LAST][currStateIndx-1]; //not transitioning, starting fr is same as before
         }
         else
         {
            dC[CURR][currStateIndx] =  dC[LAST][currStateIndx]+ LOG_0_9 + lkld[cl];
            tC[CURR][currStateIndx] = tC[LAST][currStateIndx]; //not transitioning, starting fr is same as before
         }
      }


      ////////DONE TILL HERE/////////////////////


      //Find the maximum likelihood over all the final states
      MaxClass = 0;
      MaxLkld = dC[CURR][MIN_DURATION-1];
      for(int cl=1; cl<nClass; cl++)
      {
         double currLkld = dC[CURR][cl*MIN_DURATION+MIN_DURATION-1];
         if((currLkld> MaxLkld) && U[cl])
         {
            MaxClass = cl;
            MaxLkld = currLkld;
         }
      }
      R[frame] = MaxClass; // speaker label
      F[frame] = tC[CURR][MaxClass*MIN_DURATION+MIN_DURATION-1]; // That is the backtracking ptr

      for(st=0;st <MIN_DURATION*nClass; ++st)
      {
         dC[LAST][st] = dC[CURR][st];
         tC[LAST][st] = tC[CURR][st];
      }

      /********print max 
       */
      //DEBUG
      //printf ("R = %d ", R[frame]);
      //printf ("F = %d\n", F[frame]);
      //if(frame > 400)
      //exit(0);

   }/*for(int frame = 1; frame <nSamples; frame++)*/

   *sent = new chunk;
   (*sent)->next = NULL;
   (*sent)->to = nSamples-1;
   (*sent)->cluster = R[(*sent)->to];
   (*sent)->prev = NULL;
   (*sent)->from = F[(*sent)->to];
   while((*sent)->from > 1)
   {
      (*sent)->prev = new chunk;
      (*sent)->prev->next = (*sent);
      (*sent) = (*sent)->prev;
      (*sent)->to = (*sent)->next->from - 1;
      //printf("cluster (*sent)->to = %d, we have to calc F[that]and R[that]\n",(*sent)->to); 
      if((*sent)->to < 0) (*sent)->to = 0;

      (*sent)->from = F[(*sent)->to];
      (*sent)->cluster = R[(*sent)->to];
      (*sent)->prev = NULL;
   }
   float viterbi_score =dC[CURR][R[nSamples-1]*MIN_DURATION+MIN_DURATION-1];
   printf("Segmentation Over with Viterbi score = %f\n", viterbi_score);
   fflush(stdout);
   delete []dC[LAST];delete []dC[CURR];delete []tC[LAST];delete []tC[CURR];
   delete []F; delete []R;delete []data;
   return viterbi_score;
}

void HMM::Calculate_Lkld(float *vector)
{
   for(int cl=0; cl<nClass; cl++)
   {
      if(U[cl])
         lkld[cl] = gmm[cl].Log_Likelihood(vector);
   }
}

void HMM::Train_HMM(FeatFileReader& fileReader, chunk *first)
{
   float **data;
   int count;
   for(int cl=0; cl<nClass; cl++)
   {
      if(U[cl])
      {
         Get_Data(fileReader, &data, first, cl, cl, &count);
         if(count > gmm[cl].M)
         {
            gmm[cl].lkld = gmm[cl].Train(data, count);
            free(*data); free(data);
         }/*if(count[cl] > gmm[cl].M)*/
         else
         {
            U[cl] = 0;
         }

      }/*if(U[cl])*/
   }/*for(int cl=0; cl<nCLass; cl++)*/
}

void HMM::Get_Data(FeatFileReader& fileReader ,float ***pData, chunk *first, int cd1, int cd2, int *pCount)
{
   *pCount = 0;
   for(chunk *ch = first; ch!=NULL; ch= ch->next)
   {
      if(ch->cluster == cd1 || ch->cluster == cd2)
      {
         (*pCount) += (ch->to-ch->from+1);
      }
   }
   *pData = (float **)malloc(sizeof(float *)*(*pCount));
   **pData =  (float *)malloc(sizeof(float)*(*pCount)*DIM);
   for(int i=1; i < *pCount; i++)
   {
      (*pData)[i] = (*pData)[i-1] + DIM;
   }
   if(*pData == NULL || **pData == NULL)
   {
      printf("Memory could not be allocated\n");
      exit(0);
   }
   int index = 0, fr_index=0;

   for(chunk *ch = first; ch!=NULL; ch= ch->next)
   {
      if(ch->cluster == cd1 || ch->cluster == cd2)
      {
         for(int smp=ch->from; smp <= ch->to; smp++)
         {
            fileReader.readfeat(fr_index,(*pData)[index]);
            index++;
            fr_index++;
         }
      }
      else
      {
         fr_index += (ch->to - ch->from +1);
      }
   }
}

void HMM::CMS(FeatFileReader& fileReader, chunk *first)
{
   int *count = new int[nClass];
   float *mfcc = new float[DIM];
   Initialize2DFloat(Clusters_Means, nClass, DIM);
   Initialize1DInt(count, nClass);
   for(chunk *ch = first; ch!=NULL; ch= ch->next)
   {
      count[ch->cluster] += (ch->to - ch->from +1);
      for(int smp=ch->from; smp <= ch->to; smp++)
      {
         fileReader.readfeat(smp, mfcc);
         //fread(mfcc, sizeof(float), DIM, fin);
         for(int d=0;d<DIM;d++)
         {
            Clusters_Means[ch->cluster][d] += mfcc[d];
         }
      }/*for(int smp=ch->from; smp < ch->to; smp++)*/
   }/*for(chunk *ch = *pFirst; ch!=NULL; ch= ch->next)*/
   for(int cl=0;cl<nClass;cl++)
   {
      if(count[cl])
      {
         for(int d=0;d<DIM;d++)
         {
            Clusters_Means[cl][d]/=(float)count[cl];
         }
      }
      else U[cl]=0;
   }
   delete []mfcc;
   delete []count;
}


int HMM::if_possible(int cd1, int cd2)
{
   for(int i=1;i<=NP[cd1][0];i++)
   {
      if(NP[cd1][i] == cd2)
         return 0;
   }
   return 1;
}

void HMM::Write(char *fname, chunk *first)
{
   FILE *fout = fopen(fname, "w");
   if(fout == NULL)
   {
      printf("Output file could not be opened\n");
      exit(0);
   }
   for(chunk *ch = first; ch!=NULL; ch=ch->next)
   {
      fprintf(fout, "from %d to %d: cluster %d\n", ch->from, ch->to, ch->cluster);
   }
   fclose(fout);
}

//This function write the RTTM file to the designated output file
// First the indices are mapped to the actual frame numbers
void HMM::WriteRttm(const char *fname, chunk *first, FeatFileReader& fileReader, int tot_num_frames, const char *mtgid)
{
   int i, fr; 
   int cluster_start, cluster_end, cluster_id;


   FILE *fout = fopen(fname, "w");
   if(fout == NULL)
   {
      printf("Output file %s could not be opened\n", fname);
      exit(0);
   }
   // map the frames back to original indices..

   int *allframes = new int[tot_num_frames];
   for (i=0; i<tot_num_frames;++i)
      allframes[i] = -1; 

   for( chunk *ch = first; ch!=NULL; ch=ch->next)
   {
      for(i=ch->from;i<=ch->to;i++) {
         fr = fileReader.findFrameIndex(i);//get mapping to original fr number
         if (fr>0) allframes[fr] = ch->cluster;
      }
   }

   int s_rate = fileReader.getSampleRate();

   //Print the last cluster
   double htk_rate = 10000000.0;

   cluster_start = 0; 
   cluster_end = 0;
   cluster_id = allframes[0];
   printf("The htk file sample rate %d \n", fileReader.getSampleRate());
   for(fr=1; fr<tot_num_frames; ++fr)
   {
      if (allframes[fr] == cluster_id)
      {
         cluster_end = fr;
      }
      else
      {
         //print cluster_start, cluster_end, cluster_id
         if(cluster_id != -1) //Not non-speech
         {
            float start_time = double(cluster_start ) * double(s_rate)/htk_rate;
            float dur_time = double(cluster_end -cluster_start+1) * double(s_rate)/htk_rate;
            fprintf(fout, "SPEAKER %s 1 %.2f %.2f <NA> <NA> %s_spkr_%d <NA>\n",mtgid, start_time, dur_time, mtgid, cluster_id);
         }
         cluster_start = fr; 
         cluster_end = fr;
         cluster_id = allframes[fr];
      }


   }
   //
   //Print the last cluster
   if(cluster_id != -1) //Not non-speech
   {
      float start_time = double(cluster_start ) * double(s_rate)/htk_rate;
      float dur_time = double(cluster_end -cluster_start+1) * double(s_rate)/htk_rate;
      fprintf(fout, "SPEAKER %s 1 %.2f %.2f <NA> <NA> %s_spkr_%d <NA>\n",
                    mtgid, start_time, dur_time, mtgid, cluster_id);
   }
   fclose(fout);
   delete [] allframes;

}

void HMM::WriteLkld(char *fname, char *name_dat)
{
   FILE *fin = fopen(name_dat, "rb");
   FILE *fout = fopen(fname, "w");
   float *vector = new float[DIM];
   int nSamples;
   fprintf(fout, "float32\n");
   int n_usable=0;
   for(int cl=0;cl<nClass;cl++)
   {
      if(U[cl])
         n_usable++;
   }
   fread(&nSamples, sizeof(int), 1, fin);
   fprintf(fout, "%d %d\n", n_usable, nSamples);
   for(int i=0;i<nSamples; i++)
   {
      fread(vector, sizeof(float), DIM, fin);
      Calculate_Lkld(vector);
      for(int cl=0; cl<nClass;cl++)
      {
         if(U[cl])
            fwrite(&lkld[cl], sizeof(float), 1, fout);
      }
   }
}

void HMM::WriteParameters(char *fname)
{
   FILE *fout = fopen(fname, "w");
   int n_usable = 0;
   for(int cl=0;cl<nClass;cl++)
   {
      if(U[cl])
         n_usable++;
   }
   fprintf(fout, "%d\n", n_usable);
   for(int cl=0;cl<nClass;cl++)
   {
      if(U[cl])
      {
         for(int m=0;m<gmm[cl].M;m++)
         {
            fprintf(fout, "%f\n", gmm[cl].Weight[m]);
            const float *gMean = gmm[cl].gArray[m].GetMean();
            for(int d=0;d<DIM;d++)
            {
               fprintf(fout, "%e\n", gMean[d]);
               fprintf(fout, "%e\n", 1.00);
            }
         }
      }
   }
   fclose(fout);
}
}

