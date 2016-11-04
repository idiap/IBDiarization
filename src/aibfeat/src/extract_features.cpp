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


#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

#include "global.h"
#include "hmmreadwrite.h"
#include "scpread.h"
#include "featconfig.h"
#include "extract_features.h"

#define PRECISION 1.0e7

double find_log_likelihood(float *frame_vec,double *meanvec,float* logVars, double *var_vec, bool *isVar, int dim);

void compute_features(ConfigVars &configVars, vector<vector< float > > &postmat)
{

  vector<FeatStream>::iterator featIter;

   __uint32_t n_frames;
   __uint32_t s_rate;
   __uint16_t  byte_per_frame;
   __uint16_t  tc;
   __uint16_t  frame_dim;

   float **frame_vec; //vector to hold one frame of speech
   double *mean_vec;  //vector to hold likelihoods of one frame of speech
   double *var_vec;  //vector to hold likelihoods of one frame of speech
   double *means;     //local means in each segment 
   double *wts;     //local weights in each segment 
   bool  *is_var;  // 
   unsigned int mix_index;
   unsigned int temp_offset;

   int i, j; //loop indices

   // Reading the scp file
   //parameters from the config : configVars.scpfile
   int *seg_ptr; //Segmentation info array 
   int *pur_seg_ptr; //speech segments with max duration as set in configVars.max_dur
   int num_segs;  //number of segments in the scp file
   unsigned int num_gmm_comps; // Number of GMM components
   unsigned int seg_len = 0, index_times_dim = 0;

   // First create segments of maximum max_dur long 
   //Required Parameters configVars.m_maxDur
   // -- read the scp file
   cout << " Reading and Processing the scp file \n";
   read_scp_file(configVars.scpfile.c_str(), &seg_ptr, &num_segs);
   unsigned int last_speech_frame = seg_ptr[2*num_segs-1];

   //string object to keep the file names
   stringstream ss(stringstream::in|stringstream::out); 

   int split;
   //Opening the output seg file file 

   //Just making sure there is enough memory: will release it soon
   pur_seg_ptr= (int*)malloc(2*(num_segs+last_speech_frame/configVars.m_maxDur)*sizeof(int)); //Do not worry about the inside huge value


   int num_of_seg_purity =0, li, lr;
   __uint32_t num_vec = 0;
   for(i=0; i < num_segs; ++i)
   {

      li = 2*i;
      lr = li + 1;
      seg_len = seg_ptr[lr] - seg_ptr[li]+1;
      num_vec += seg_len;
      split = floor(seg_len/configVars.m_maxDur);
      if (split >= 1)
      {
         long beg_purity=seg_ptr[2*i], end_purity=seg_ptr[2*i]+configVars.m_maxDur-1;
         for(j=0; j < split; j++)
         {
            pur_seg_ptr[2*num_of_seg_purity] = beg_purity;
            pur_seg_ptr[2*num_of_seg_purity+1] = end_purity;
            num_of_seg_purity++;
            beg_purity = end_purity+1;
            end_purity = end_purity+configVars.m_maxDur;
         }
         if(beg_purity <=  seg_ptr[2*i+1] )
         {
            pur_seg_ptr[2*num_of_seg_purity] = beg_purity;
            pur_seg_ptr[2*num_of_seg_purity+1] = seg_ptr[2*i+1] ;
            num_of_seg_purity++;
         }
      }
      else
      {
         pur_seg_ptr[2*num_of_seg_purity] = seg_ptr[2*i] ;
         pur_seg_ptr[2*num_of_seg_purity+1] = seg_ptr[2*i+1] ;
         num_of_seg_purity++;
      }
   }

   cout << "Number of vectors = " << num_vec << endl;
   //Copy pur_seg_pointer to seg_ptr
   num_segs = num_of_seg_purity;
   free(seg_ptr); seg_ptr = NULL;

   //Making both of them the same
   seg_ptr = pur_seg_ptr;

   //Opening an scp file 
   ss << configVars.m_tmpDir << "/"<<configVars.id<<".scp";
   FILE *segfp = fopen(ss.str().c_str(),"w");
   if(segfp == NULL )
   {
      cerr<<"Unable to write segfp file \""<<ss.str().c_str()<<"\""<<endl;
      exit(0);
   }
   ss.str(""); //Flush that buffer  for other things

   for(i =0; i <num_segs; ++i)
   {
      int start = seg_ptr[2*i];
      int end = seg_ptr[2*i+1];
      fprintf(segfp,"%s_%d_%d=%s.scp[%d,%d]\n",configVars.id.c_str(),start,end,configVars.id.c_str(),start,end);
   }
   fclose(segfp);

   /////////////////////////////////////////////////////
   unsigned int diar_feature_index=0;
   unsigned int f_index, d_index; //feature(time) and dim index

   
   for (featIter = configVars.m_feats.begin(); featIter < configVars.m_feats.end(); ++featIter)
   {

      const char *curr_feat_file = featIter->m_filename.c_str();
      FeatType currFeatType = featIter->m_feattype;

      cout << "Feat file "<<curr_feat_file<< " of type "<< currFeatType <<endl;
      //Open and read the file 
      FILE *featfp = fopen(curr_feat_file,"rb");

      printf("Reading feature file %s\n",curr_feat_file);
      if (featfp == NULL)
      {
         fprintf(stderr,"Could not open feature file for reading <%s>\n",curr_feat_file);
         exit(1);
      }
      //Reading HTK Header
      fread(&n_frames,SIZE_LONG,1,featfp); //Number of frames
      fread(&s_rate,SIZE_LONG,1,featfp); //Sampling Rate
      fread(&byte_per_frame,SIZE_SHORT,1,featfp); //Number of byte_per_frame per frame
      fread(&tc,SIZE_SHORT,1,featfp); //tc ??

      n_frames = SwapByteOrderOfLong(&n_frames);
      byte_per_frame = SwapByteOrderOfShort(&byte_per_frame);
      frame_dim = byte_per_frame / SIZE_LONG;
      cout << "frame dim: " << frame_dim << endl;

      frame_vec = NULL;
      cout << "num_vec = " << num_vec << endl;      
      frame_vec = (float**)malloc((n_frames)* sizeof(float*));
      if (frame_vec == NULL) {
          cout << "Unable to allocate large amounts of memory." << endl;
          exit(2);
      }
      unsigned int ii;
      cout << "Attempting memory allocation." << endl;
      fflush(stdout);
      for(ii = 0; ii < n_frames; ii++) {
          frame_vec[ii] = NULL;
          frame_vec[ii] = (float *) malloc(frame_dim * sizeof(float));
          if (frame_vec[ii] == NULL) {
              cout << "Unable to allocate frame." << endl;
              exit(2);
          }
      }
      cout << "Memory successully allocated." << endl;
      fflush(stdout);
      mean_vec  = (double*)malloc(frame_dim * sizeof(double));
      var_vec   = (double*)malloc(frame_dim * sizeof(double));
      is_var    = (bool*)malloc(frame_dim * sizeof(bool));
      double *prec_vec = (double *) malloc(frame_dim * sizeof(double));

      //Part I - Find feats with low Variance / Also compute local mean vectors
      means = (double *)malloc(num_segs*frame_dim*sizeof(double));
      wts = (double *)malloc(num_segs*sizeof(double));
      memset(means,0,num_segs*frame_dim*sizeof(double));
      memset(wts,0,num_segs*sizeof(double));

      //Initializing the accumulators
      for(i=0; i<frame_dim; ++i)
      {
         mean_vec[i]=0.0;
         var_vec[i]=0.0;
      }
      if (fseek(featfp,2*(SIZE_LONG+SIZE_SHORT),SEEK_SET) != 0) 
          std::cout << "ERROR!!!!" << std::endl;

      float frame_val = 0;
      int seg_index = 0;
      unsigned int first_index = seg_ptr[0]-1;
      unsigned int last_index  = seg_ptr[1]-1;
      int num_sp_frames = 0;
      //Accumulating the mean and variance statistics
      unsigned int idx = 0;
      cout << "Reading file ... " << endl;
      fflush(stdout);
      for(f_index=0; f_index < n_frames; ++f_index)
      {
         if(fread(frame_vec[f_index],sizeof(float),frame_dim,featfp) == 0) {
             cout << "Reached end of file." << endl;
             break;
         }
         unsigned char uiter, *yiter = NULL;

         for(d_index=0; d_index < frame_dim; ++d_index) {
            yiter = (unsigned char *) &frame_vec[f_index][d_index];
            uiter = yiter[0];
            yiter[0] = yiter[3];
            yiter[3] = uiter;
            uiter = yiter[1];
            yiter[1] = yiter[2];
            yiter[2] = uiter;
         }
         if (f_index < first_index || f_index > last_index) 
           continue;

         //Only considering the speech frames
            for(d_index=0; d_index < frame_dim;  ++d_index)
            {
               frame_val = frame_vec[f_index][d_index];
               if(isnan(frame_val)) {
                   cout <<"Warning frame_val is nan for idx = " << f_index << " and d_index == " << d_index << " f_index = " << f_index << endl;
                   exit(1);
               }
               mean_vec[d_index] += frame_val;
               var_vec[d_index] += (frame_val * frame_val);
               means[index_times_dim + d_index] += frame_val;
            }
            if(f_index == last_index)
            {
               seg_len = last_index-first_index +1;
               for(d_index=0; d_index<frame_dim; ++d_index)
                  means[index_times_dim + d_index] /= seg_len;
               wts[seg_index] = seg_len;
               num_sp_frames += wts[seg_index];

               ++seg_index;
               index_times_dim += frame_dim;
               if(seg_index < num_segs) {
                  first_index = seg_ptr[2*seg_index]-1;
                  last_index = seg_ptr[2*seg_index+1]-1; 
                  if(last_index >= n_frames) last_index = n_frames-1; 
                  fflush(stdout);
               }
               else {
                   cout << "all segments read " << endl;
                   fflush(stdout);
                   break;
               }
            }
      }

      //num_segs=seg_index; //In case the number of segments is less

      cout << "Number of segments in the file = " << num_segs << endl;
      int nfeats = 0;

      // performing feature selection and 
      //variance smoothing 
      for(i=0; i<num_segs; ++i)
         wts[i] /=  num_sp_frames;

      float smooth_val = (currFeatType == TDOA)?60.0:25.0; // The smooth value is hard-coded for TDOA features
      for(i=0; i<frame_dim; ++i)
      {
         mean_vec[i]/=num_sp_frames;
         var_vec[i]/=num_sp_frames;
         var_vec[i]-=(mean_vec[i]*mean_vec[i]);
         is_var[i] = (var_vec[i]>VAR_TH);
         if(is_var[i])
            ++nfeats;
         
         //performing smoothing 
         var_vec[i] /= smooth_val;
         if(var_vec[i] < 10E-6) 
             var_vec[i] = 10E-6;
         prec_vec[i] = 1.0 / var_vec[i];
      }


      FILE *outfp = NULL;
      num_gmm_comps = seg_index;

      // printing the HMM to a file
      {
         double *iter; 

         ss << configVars.m_tmpDir << "/"<<configVars.id<<".hmm."<< diar_feature_index;
         outfp=fopen(ss.str().c_str(),"w");
         ss.str(""); //Flush that buffer  for other things
         fprintf(outfp,"~o\n<STREAMINFO> 1 %d\n<VECSIZE> %d<NULLD><MFCC><DIAGC>\n",frame_dim, frame_dim);
         fprintf(outfp,"~h \"spkr\"\n<BEGINHMM>\n<NUMSTATES> 3\n<STATE> 2\n<NUMMIXES> ");
         fprintf(outfp,"%d \n",num_segs); 
         iter = var_vec; 

         for (i=0; i < num_segs; i++)
         {
            fprintf(outfp,"<MIXTURE> %d %f\n<MEAN> %d\n ",i+1,wts[i],frame_dim);

            iter =  means + frame_dim*i;
            for(j=0; j < frame_dim; ++j,++iter)
               fprintf(outfp,"%f ",*iter);

            fprintf(outfp,"\n<VARIANCE> %d\n",frame_dim);

            iter = var_vec; 
            for(j=0; j < frame_dim; ++j,++iter)
               fprintf(outfp,"%f ",*iter);

            //for(j=0; j < frame_dim; ++j,++iter)
            //fprintf(stderr,"%f ",*iter/smooth_val);
            //fprintf(stderr,"\n");

            fprintf(outfp,"\n");
            fprintf(outfp,"<GCONST> 1\n");
         }


         fprintf(outfp,"<TRANSP> 3\n0.000000e+00 1.000000e+00 0.000000e+00\n"
                       "0.000000e+00 9.953219e-01 4.678072e-03\n0.000000e+00"
                       "0.000000e+00 0.000000e+00\n<ENDHMM>");
         fclose (outfp);
      }


      //Part II - estimate the GMM 
      //Required Parameters configVars.m_maxDur
      //global var and local means computed in the previous step
      // and smoothing also in the next step

      // Part - III 
      //features to Posteriors 
      if (fseek(featfp,2*(SIZE_LONG+SIZE_SHORT),SEEK_SET) != 0 ) std::cout << "ERROR!!!!!" << std::endl;
      float *lkld = (float*) malloc(num_gmm_comps * sizeof(float));
      float *logwts = (float*) malloc(num_gmm_comps * sizeof(float));
      float *logVars = (float*) malloc(frame_dim * sizeof(float));
      double *postrs  = (double*) malloc(num_gmm_comps * sizeof(double));
      double l2_x;
      float *segfeat = (float*) malloc(num_gmm_comps * sizeof(float));

      memset(lkld,0,num_gmm_comps*sizeof(float));
      memset(postrs,0,num_gmm_comps*sizeof(double));
      memset(segfeat,0,num_gmm_comps*sizeof(float));
      temp_offset = 0;
      for (mix_index = 0; mix_index < num_gmm_comps; ++mix_index) {
         logwts[mix_index] = log(wts[mix_index]);
         
         l2_x = 0.;
         for (j=0; j<frame_dim; j++) {
             l2_x += (means[temp_offset + j] * means[temp_offset + j] * prec_vec[j]);
         }
         logwts[mix_index] = logwts[mix_index] - 0.5 * l2_x;
         if(isnan(logwts[mix_index])) {
             cout << "Logwts is nan for mix_index " << mix_index << " and featuretype " << currFeatType << endl;
             cout << " logwts = " << logwts[mix_index] << " and wt = " << wts[mix_index] << " l2_x = " << l2_x << endl;
             exit(0);
         }
         temp_offset += frame_dim;         
      }
      temp_offset = 0;
      for(mix_index = 0; mix_index < num_gmm_comps; ++mix_index) {
          for(j = 0; j < frame_dim; j++) {
              means[temp_offset + j] *= prec_vec[j];
          }
          temp_offset += frame_dim;
      }

      //posterior output file
      ss << configVars.m_tmpDir << "/"<<configVars.id<<".post."<< diar_feature_index; 
      outfp=fopen(ss.str().c_str(),"wb");
      ss.str(""); //Flush that buffer  for other things
      if(outfp == NULL)
      {
         cout << "Cannot open output file "<<ss.str().c_str()<<endl; 
         exit(0);
      }


      //Segment wise accumulated posteriors 
      ss << configVars.m_tmpDir << "/"<<configVars.id<<".feat."<< diar_feature_index;
      FILE *segfp=fopen(ss.str().c_str(),"w");
      ss.str(""); //Flush that buffer  for other things
      if(outfp == NULL)
      {
         cout << "Cannot open output file "<<ss.str().c_str()<<endl; 
         exit(0);
      }


      __uint32_t n_frames_rev = SwapByteOrderOfLong(&n_frames);
      __uint16_t byte_per_frame_rev = SwapByteOrderOfShort(&byte_per_frame);

      tc = 9; //USER Datatype in HTK
      tc = SwapByteOrderOfShort(&tc);
      byte_per_frame_rev = num_gmm_comps * SIZE_LONG;
      byte_per_frame_rev = SwapByteOrderOfShort(&byte_per_frame_rev);
      fwrite(&n_frames_rev,SIZE_LONG,1,outfp); //Number of frames
      fwrite(&s_rate,SIZE_LONG,1,outfp); //Sampling Rate
      fwrite(&byte_per_frame_rev,SIZE_SHORT,1,outfp); //Number of byte_per_frame per frame
      fwrite(&tc,SIZE_SHORT,1,outfp); //tc ??

      seg_index = 0;
      first_index = seg_ptr[0];
      last_index  = seg_ptr[1];

      // II time looping through the file
      float *post_zero = (float *) malloc(sizeof(float)*num_gmm_comps);
      for(mix_index = 0; mix_index < num_gmm_comps; mix_index++) {
          post_zero[mix_index] = 0.;
          post_zero[mix_index] = SwapByteOrderOfFloat(&post_zero[mix_index]);
      }

      cout << "num_vec = " << num_vec << "  idx = " << idx << endl;
      for(f_index = 0; f_index< n_frames ; ++f_index)          
      {

         if (f_index < first_index || f_index > last_index) {
             fwrite(post_zero, 4, num_gmm_comps, outfp);
             continue;
         }
         float maxLkld = 0;

         //Finding the log likelihood of each frame (is_var is supposed to take care of eliminated dimensions)
         int max_index = 0; 
         double temp_term2;
         
         for(mix_index = 0; mix_index < num_gmm_comps; ++mix_index)
         {
            double *meanvec = means+ mix_index*frame_dim ;
            temp_term2 = logwts[mix_index];
            if(isnan(temp_term2)) {
                cout << "logwt is nan for idx " << f_index << " ,mix_index " << mix_index << " feature " << currFeatType << endl;
                exit(1);
            }
            for(j = 0; j < frame_dim; j++) {
                temp_term2 += frame_vec[f_index][j] * meanvec[j]; 
            }
            lkld[mix_index] = temp_term2;
            if(isnan(lkld[mix_index])) {
                cout << "lkld is nan for idx " << f_index << " ,mix_index " << mix_index << " feature " << currFeatType << endl;
                exit(1);
            }

            if ((lkld[mix_index] > maxLkld || mix_index == 0 )) {maxLkld = lkld[mix_index]; max_index = mix_index;}
         }
         idx++;

         double sumLkld = 0;
         //Subtracting the max likelihood and exp (also finding the normalizing term)
         for (mix_index = 0; mix_index < num_gmm_comps; ++mix_index)
         {
            lkld[mix_index] -= maxLkld;
            postrs[mix_index] = exp(double(lkld[mix_index]));
            sumLkld += postrs[mix_index];
         }

         //Normalizing and finding the max (threshold)
         for (mix_index = 0; mix_index < num_gmm_comps; ++mix_index)
         {
            postrs[mix_index]/=sumLkld;
            float postrs_rev = postrs[mix_index]; 
            postrs_rev = SwapByteOrderOfFloat(&postrs_rev);
            fwrite(&postrs_rev,sizeof(float),1,outfp);
         }
         if( (f_index>= first_index) && (f_index <= last_index))
         {
            if(currFeatType == TDOA) //TDOA have more peaky distributions
            {
               for(d_index=0; d_index < num_gmm_comps;  ++d_index)
                  segfeat[d_index] += postrs[d_index];
            }
            else  
            {
               if(max_index > -1)
                  segfeat[max_index] += 1.0 ; 
            }
            if(f_index == last_index)
            {
               for(d_index=0; d_index < num_gmm_comps; ++d_index)
               {
                  fprintf(segfp,"%g ",segfeat[d_index]);
                  segfeat[d_index] = 0;
               }
               fprintf(segfp,"\n");
               ++seg_index;
               if(seg_index < num_segs) {
                  first_index = seg_ptr[2*seg_index];
                  last_index = seg_ptr[2*seg_index+1];
                  if(last_index >= n_frames) last_index=n_frames-1;
               }
            }
         }
      }// Finish looping through the frames for II time
      fclose (outfp);
      fclose (segfp);
      fclose (featfp);

      free(lkld);
      free(logwts); 
      free(logVars); 
      free(postrs); 
      free(segfeat); 

      //seg feat
      ++diar_feature_index;
      for(ii = 0; ii < last_speech_frame && ii << n_frames; ii++) {
          if(frame_vec[ii] != NULL)
              free(frame_vec[ii]);          
      }
      seg_index = 0;
      index_times_dim = 0;
      free (frame_vec);
      free (mean_vec);
      free (var_vec);
      free (prec_vec);
      free (is_var);
      free (means);
      free (wts);

   } // Finished Looping through all features

   //Allocating two arrays of file pointers for posterior and feat files
   int numFeats = configVars.m_feats.size();
   FILE **inpostfp = new FILE*[numFeats]; //input posterior file pointers for each feature
   FILE **infeatfp = new FILE*[numFeats]; //input acc.postr.file pointers for each feature
   FILE *outpostfp;
   FILE *outfeatfp;

   int featIdx = 0;

   ss.str(""); //Flush that buffer  for other things
   //Opening the output posterior file
   ss << configVars.m_tmpDir << "/"<<configVars.id<<".post.all";
   outpostfp = fopen(ss.str().c_str(),"wb");
   if(outpostfp ==NULL)
   {
      cerr<<"Unable to open output file \""<<ss.str().c_str()<<"\" to accumulate all posteriors"<<endl;
      exit(0);
   }
   ss.str(""); //Flush that buffer  for other things

   ss << configVars.m_tmpDir << "/"<<configVars.id<<".feat.all";
   outfeatfp = fopen(ss.str().c_str(),"w");
   if(outfeatfp ==NULL)
   {
      cerr<<"Unable to open output file \""<<ss.str().c_str()<<"\" to accumulate all features"<<endl;
      exit(0);
   }
   ss.str(""); //Flush that buffer  for other things

   int seg_index = 0;
   frame_dim = 0; 
   //combine features  Open an array of file pointers for each posterior
   //and accumulated feature files (input to clustering)
   for (featIter = configVars.m_feats.begin(),featIdx = 0; featIter < configVars.m_feats.end(); ++featIdx, ++featIter)
   {

      //Opening posterior file
      ss << configVars.m_tmpDir << "/"<<configVars.id<<".post."<< featIdx;
      inpostfp[featIdx]=fopen(ss.str().c_str(),"rb");
      if(inpostfp[featIdx] == NULL )
      {
         cerr<<"Unable to open input feat file \""<<ss.str().c_str()<<"\" of index "<<featIdx<<endl;
         exit(0);
      }
      ss.str(""); //Flush that buffer  for other thing
      if (featIdx == 0) {
         fread(&n_frames,SIZE_LONG,1,inpostfp[featIdx]); //Number of frames
         fread(&s_rate,SIZE_LONG,1,inpostfp[featIdx]); //Sampling Rate
         fread(&byte_per_frame,SIZE_SHORT,1,inpostfp[featIdx]); //Number of byte_per_frame per frame
         fread(&tc,SIZE_SHORT,1,inpostfp[featIdx]); //tc ??

         fwrite(&n_frames,SIZE_LONG,1,outpostfp); //Number of frames
         fwrite(&s_rate,SIZE_LONG,1,outpostfp); //Sampling Rate
         fwrite(&byte_per_frame,SIZE_SHORT,1,outpostfp); //Number of byte_per_frame per frame
         fwrite(&tc,SIZE_SHORT,1,outpostfp); //tc ??

         n_frames = SwapByteOrderOfLong(&n_frames);
         byte_per_frame = SwapByteOrderOfShort(&byte_per_frame);
         frame_dim = byte_per_frame / SIZE_LONG;
      }
      else{
         if ( fseek(inpostfp[featIdx],2*(SIZE_LONG+SIZE_SHORT),SEEK_SET) != 0 ) std::cout << "ERROR!!!!!" << std::endl;
      }


      //Opening accumulated posteriors (input to clustering)
      ss << configVars.m_tmpDir << "/"<<configVars.id<<".feat."<< featIdx;
      infeatfp[featIdx]=fopen(ss.str().c_str(),"r");
      if(infeatfp[featIdx] == NULL )
      {
         cerr<<"Unable to open input feat file \""<<ss.str().c_str()<<"\" of index "<<featIdx<<endl;
         exit(0);
      }
      ss.str(""); //Flush that buffer  for other things
   }// Finished Looping through all features

   //Allocate the frame memory
   float *in_frame = new float [frame_dim]; 
   float *out_frame = new float [frame_dim];

   //Step I - Go through each frame of posteriors and sum them
   for(f_index=0; f_index<n_frames; ++f_index)
   {
      //initialize the accumulator to zero
      for(d_index = 0; d_index < frame_dim; ++d_index)
         out_frame[d_index] = 0;

      float *iter = NULL; //Temp pointer for inp/out byte reverse
      for (featIter = configVars.m_feats.begin(),featIdx = 0; featIter < configVars.m_feats.end(); ++featIdx, ++featIter)
      {
         fread(in_frame,SIZE_LONG,frame_dim,inpostfp[featIdx]);
         float feat_wt = featIter->m_wt;

         //weighted addition
         for(d_index=0,iter=in_frame; d_index < frame_dim;++iter, ++d_index)
            out_frame[d_index] += feat_wt*SwapByteOrderOfFloat(iter);
      }

      for(d_index=0,iter=out_frame; d_index < frame_dim;++iter, ++d_index)
         *iter = SwapByteOrderOfFloat(iter);

      fwrite(out_frame,SIZE_LONG,frame_dim,outpostfp);
   }

   double* out_acc_frame = new double[frame_dim];

   for (seg_index=0; seg_index < num_segs; ++seg_index)
   {
      //initialize the accumulator to zero
      for(d_index = 0; d_index < frame_dim; ++d_index)
         out_acc_frame[d_index] = 0;

      //Read one segment from all the files and accumulate
      for (featIter = configVars.m_feats.begin(),featIdx = 0; featIter < configVars.m_feats.end(); ++featIdx, ++featIter)
      {
         float feat_wt = featIter->m_wt;
         float *iter = in_frame;
         for(d_index=0; d_index < frame_dim; ++d_index)
         {
            fscanf(infeatfp[featIdx],"%g", iter);

            //weighted addition
            out_acc_frame[d_index] += feat_wt*(*iter);

         }
      }
      //writing to the output
      
      for(d_index = 0; d_index < frame_dim; ++d_index)
	{ 
	  fprintf(outfeatfp,"%g ",out_acc_frame[d_index]); 
	}
      
      postmat.push_back(vector<float> (out_acc_frame,out_acc_frame + frame_dim) );
      
      fprintf(outfeatfp,"\n");
   }

   for (featIdx = 0; featIdx < numFeats ; ++featIdx)
   {
      fclose(inpostfp[featIdx]);
      fclose(infeatfp[featIdx]);
   }
   fclose (outpostfp);
   fclose (outfeatfp);

   delete [] in_frame;
   delete [] out_frame;
   delete [] out_acc_frame;

   delete [] inpostfp;
   delete [] infeatfp;
   free (pur_seg_ptr);  
}

double find_log_likelihood(float *frame_vec,double *meanvec,float* logVars, double *prec_vec, bool *isVar, int dim)
{
   double likelihood = 0;
   double *mean, *var;
   float *frame;
   int i;
   for (i = 0, mean = meanvec, var=prec_vec, frame=frame_vec;  i< dim ; ++i,++mean,++var,++frame)
   {
      if(isVar[i]){
         double lkld =  (*frame - *mean); 
         lkld = lkld*lkld*(*var);
         likelihood += lkld;
      }
   }
   likelihood = -likelihood*0.5 ;
   return likelihood;
}

