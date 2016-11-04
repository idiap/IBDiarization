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
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <sstream>


using namespace std;
#include "tools.h"

void readmatrixfromfile(const char * infilename, vector <vector <double> > *MM)
{
   string line;  
   ifstream myfile (infilename); 
 
   unsigned int indx=0;
   if (myfile.is_open()) {
       while (myfile.good()) {
           getline (myfile,line);
           
           std::istringstream stm(line);
           istream_iterator<double> begin(stm), end ;
           std::vector<double> values(begin, end);
           
           if (values.size() > 0) {
               copy(values.begin(),values.end(),(*MM)[indx].begin());
               indx++;
           }
       }
       myfile.close();
       cout << "\n\n";
   }
}


void readmatrixfromfile_and_resize(const char * infilename, vector <vector <double> > *MM)
{
   string line;  
   ifstream myfile (infilename); 
 
   unsigned int indx=0;
   if (myfile.is_open()) {
       while (myfile.good()) {
	   getline (myfile,line);
	   
	   std::istringstream stm(line);
	   istream_iterator<double> begin(stm), end ;
	   std::vector<double> values( begin, end );
	   
	   if (values.size() > 0) {
	       if (indx>0) { 
               (*MM).resize( (*MM).size()+1 );   
           }
	       
	       (*MM)[indx].resize( values.size() );
	       
	       copy(values.begin(),values.end(),(*MM)[indx].begin());
	       indx++;
	     }
       }
       myfile.close();
   }
}


void printvector_on_file(const char* outfilename, vector <int>* vout)
{
   ofstream myfile (outfilename);
   
   ostream_iterator<int> out_it (myfile," ");
   copy((*vout).begin(),(*vout).end(),out_it);
}
