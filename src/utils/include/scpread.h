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
#include <algorithm>
#include <functional>
#include <vector>
#include <set>

#define LINE_LENGTH 1024

using namespace std;

//This is for sorting the vectors 
template <class T>
class isvec_less_than {                            // A function object
   public:
      bool operator()( const T &x, const T &y ) { return x[0] < y[0];}
};


int read_scp_file(const char *in_filename, int **out_seginfo, int *num_segs);
bool is_seg_lessthan(vector<int> , vector<int> );


