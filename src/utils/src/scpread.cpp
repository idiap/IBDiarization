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


#include "scpread.h"
/* function read_scp_file
 * This function reads an scp file and 
 * returns the sorted seginfo in the 
 * float array (adjacent ones r the pairs)
 */

int read_scp_file(const char *in_filename, int **out_seginfo, int *num_segs)
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
      //printf("SCP READ %d %d\n",ivel[0], ivel[1]);
   }

   long last_speech_frame = ivel[1];

   vector<bool> is_speech_frame(last_speech_frame+1,false);
   for (unsigned int i =0; i <pos.size(); ++i)
   {
      ivel = pos[i];
      for  (int j=ivel[0]; j <ivel[1]; ++j)
         is_speech_frame[j] = true; 
   }

   pos.clear();
   bool sp_flag = false; 
   start = end = -1;
   for(long fr = 1; fr<=last_speech_frame; ++fr)
   {
      if(!sp_flag && is_speech_frame[fr])
         start = fr;
      else if ((sp_flag && !is_speech_frame[fr]) || (fr >= last_speech_frame))
      {
         end = fr;
         ivel.clear();
         ivel.push_back(start);
         ivel.push_back(end);
         pos.push_back(ivel);
         //printf(" speech seg = %d %d\n",ivel[0], ivel[1]);
      }
      sp_flag = is_speech_frame[fr];
   }

   *out_seginfo = (int*) malloc(2*pos.size()*sizeof(int));
   for (unsigned int i =0; i <pos.size(); ++i)
   {
      (*out_seginfo)[2*i] = pos[i][0];
      (*out_seginfo)[2*i+1] = pos[i][1];
      //printf("< %7d %7d >\n",(*out_seginfo)[2*i],(*out_seginfo)[2*i+1]);
   }

   //printf("\n Sorted List \n");
   //for (int i =0; i <pos.size(); ++i)
   //printf("%d %d\n",pos[i][0], pos[i][1]);

   *num_segs = pos.size();
   printf("number of segments inside = %d\n",*num_segs);
   fclose(scpfp);
   return 0;
}

   bool is_seg_lessthan(vector<int> a, vector<int> b){
      if ((a.size()>1 && b.size()>1) && (a[1] <= b[1]))
         return true;
      else 
         return false;
   }

