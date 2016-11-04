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
#include "global.h"

#define MAX_LEN 100000
#define EPS 1e-05


int read_htk_gmm(char *hmm_filename, GMM* gmm);
int write_htk_gmm(char *hmm_filename, GMM* gmm);


int read_htk_gmm(char *hmm_filename, GMM* gmm)
{
   //This function reads a gmm from a htk hmm file
   //The input should be a one state hmm in ascii
   //format
   //params
   //   hmm_filename -> input hmm filename
   //   gmm          -> gmm data structure

   int vec_size;
   int num_states =0;
   int num_models;
   int num_mixes;

   FILE *hmmfp = fopen(hmm_filename, "r");
   if(hmmfp == NULL)
   {
      fprintf(stderr,"FILE OPEN ERROR\n");
      exit(0);
   }

   char *line = (char*)malloc(MAX_LEN *sizeof(char));

   /***********************/
   // reading the model

   num_models = -1;

   //--------------------------------------------------
   //-- Calculating the #models, dim, and #mixtures  --
   //--------------------------------------------------
   fgets(line,MAX_LEN,hmmfp);
   while (!feof(hmmfp)) {
      if (strstr(line,"~h")) {
         num_models++;
      } else if (strstr(line,"<NUMSTATES>")) {
         sscanf(line,"%*s %d",&num_states);
      }
      else if (strstr(line,"<VECSIZE>")) {
         sscanf(line,"%*s %d",&vec_size);
      }
      else if (strstr(line,"<NUMMIXES>")) {
         sscanf(line,"%*s %d",&num_mixes);
      }
      fgets(line,MAX_LEN,hmmfp);
   }

   if(num_models >1){
      fprintf(stderr,"ERROR: Does not Support Reading multiple models\n");
      return -1;
   }
   if(num_states !=3){
      fprintf(stderr,"ERROR: file=<%s>\n\t #states in this hmm = %d Currenly supporting only #states=1\n",hmm_filename, num_states);
      return -1;
   }

   printf("num_mixes = %d vec_size = %d \n ", num_mixes, vec_size);
   gmm->Initialize(num_mixes,vec_size);

   //----------------------------------------------//
   //--- Reading the weights means and variances --//
   //----------------------------------------------//
   int mix_index;
   float wt;
   char *tmp_ptr;
   float *tmp_flt = NULL;
   int l;
   int val;

   tmp_flt = (float*)malloc (vec_size*sizeof(float));

   mix_index = -1;

   rewind(hmmfp);
   fgets(line,MAX_LEN,hmmfp);
   while (!feof(hmmfp)) {
      if (strstr(line,"<MIXTURE>")) {
         ++mix_index;
         sscanf(line,"%*s %d %f",&val,&wt);  //Reading the mixture index
         //and corresp. weight
         gmm->Weight[mix_index] = wt;
      } 
      else if (strstr(line,"<MEAN>")) {
         if(mix_index >-1){
            fgets(line,MAX_LEN,hmmfp);
            //tmp_flt = (gmm->gArray[mix_index-1]).Mean;
            tmp_ptr = strtok(line," ");
            for (l=0;l<vec_size;l++) {
               tmp_flt[l] = atof(tmp_ptr);
               tmp_ptr = strtok(NULL," "); 
            }
            (gmm->gArray[mix_index]).SetMean(tmp_flt);
         }
      }
      else if (strstr(line,"<VARIANCE>")) {
         if(mix_index > -1){
            fgets(line,MAX_LEN,hmmfp);
            //tmp_flt = (gmm->gArray[mix_index-1]).Var;
            tmp_ptr = strtok(line," ");
            for (l=0;l<vec_size;l++) {
               tmp_flt[l] = atof(tmp_ptr);
               tmp_ptr = strtok(NULL," "); 
            }
            (gmm->gArray[mix_index]).SetVar(tmp_flt);
         }
      }
      fgets(line,MAX_LEN,hmmfp);
   }

   free (tmp_flt);
   fclose(hmmfp);
   free(line);
   return 0;
}


int write_htk_gmm(char *hmm_filename, GMM* gmm)
{

   // This function writes a given gmm as a one state
   // hmm in htk ascii format
   //
   // params
   //    hmm_filename - output filename
   //    gmm          - input gmm data structure
   //
   //    ATTENTION!! -> The program just SKIPS components with 
   //                   negative weights.. NO WARNINGS!!!!!!!!

   int i,j;
   FILE *outfp;   //Posteriors output file pointer
   int num_mix = 0;

   for (i=0; i<gmm->M; ++i)
   {
      if(gmm->Weight[i] > 0)
         ++num_mix;
   }
   const float *vec;
   outfp = fopen(hmm_filename,"w");
   if(outfp ==NULL)
   {
      fprintf(stderr,"ERROR: Cannot open output file for writing <%s>\n",hmm_filename);
      return -1;
   }
   fprintf(outfp,"~o\n<STREAMINFO> 1 %d\n",gmm->DIM);
   fprintf(outfp,"<VECSIZE> %d <NULLD><MFCC><DIAGC>\n",gmm->DIM);
   fprintf(outfp,"~h \"spkr\"\n");
   fprintf(outfp,"<BEGINHMM>\n<NUMSTATES> 3 \n<STATE> 2\n<NUMMIXES> %d\n",num_mix);

   for (i=0; i<gmm->M; ++i)
   {
      if(gmm->Weight[i] > 0)
      {
         //weight
         fprintf(outfp,"<MIXTURE> %d %e\n",i+1,gmm->Weight[i]);

         //Mean
         fprintf(outfp,"<MEAN> %d\n",gmm->DIM);
         vec = gmm->gArray[i].GetMean();
         for (j=0; j<gmm->DIM; ++j)
            fprintf(outfp,"%e ",vec[j]);
         fprintf(outfp,"\n");

         //var
         fprintf(outfp,"<VARIANCE> %d\n",gmm->DIM);
         vec = gmm->gArray[i].GetVar();
         for (j=0; j<gmm->DIM; ++j)
            fprintf(outfp,"%e ",vec[j]);
         fprintf(outfp,"\n");
      }
   }
   fprintf(outfp,"<TRANSP> 3\n");
   fprintf(outfp,"0.0 1.0 0.0\n");
   fprintf(outfp,"0.0 0.9 0.1\n");
   fprintf(outfp,"1.0 0.0 0.0\n");
   fprintf(outfp,"<ENDHMM>\n");
   fclose(outfp);
   return 0;
}
