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
#include <math.h>
#include <memory.h>
#include <functional>
#include <numeric>
#include <iostream>
#include <vector>
#include <ctime>
#include <time.h>
#include <limits>
#include <omp.h>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm> 
#include <cctype>
#include <locale>
#include <boost/lexical_cast.hpp>
#include "aglobal.h"
#include "aIB.h"
#include "fileIO.h"
#include "tools.h"
#include "various_exceptions.h"
#include "aibfunctions.h"
using namespace std;

// trim from start
static inline string &ltrim(string &s) {
        s.erase(s.begin(), find_if(s.begin(), s.end(), not1(std::ptr_fun<int, int>(isspace))));
        return s;
}

// trim from end
static inline string &rtrim(string &s) {
        s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline string &trim(string &s) {
        return ltrim(rtrim(s));
}

int main(int argc, char *argv[]) {
    FILE            *featfile = NULL ,
                    *clustfile = NULL,
                    *scpfile = NULL;
    double          beta_tvalue, 
                    nmi_tvalue; 
    int             maxclustnum;
    string          ipfilename, 
                    opfilename,
                    result_dir,
                    meeting_id;
                    

    if(argc < 3 ) {
        cout << "Usage: runaib meeting_id config_file" << endl;
        exit(0);
    }

    meeting_id = argv[1];

    ifstream ip_file(argv[2]);
    string line;
    while( getline(ip_file, line) )
    {
      if(line.length() == 0 || line.at(0) == '#')  
          continue;
      
      istringstream is_line(line);
      string key;
      if( getline(is_line, key, '=') )
      {
        string value;
        if( getline(is_line, value) )  {
            int idx = value.find('#');
            if(idx > -1) 
                value = value.substr(0, idx-1);
            
            value = rtrim(value);
            if(key == "MAX_CLUSTERS") {
                maxclustnum = lexical_cast<int>(value);
            }
            else if(key == "MODEL_SEL_PARAM") {
                nmi_tvalue = lexical_cast<double>(value);
            }
            else if(key == "AIB_BETA") {
                beta_tvalue = lexical_cast<double>(value);
            }
            else if(key == "TMP_DIR" ) {
                result_dir = string(value);
            }
            else
                continue;
        }        
      }
    }    

    opfilename = result_dir + "/" + meeting_id + ".clust.out";
    cout << "output file name is : " << opfilename << endl;
    
    ipfilename = result_dir + "/" + meeting_id + ".feat.all";    
    cout << "input flie name is: " << ipfilename << endl;

    
    functionals thisfunc;
    thisfunc.threads_num = 1;

    load_prepare_and_cluster(ipfilename.c_str(),  opfilename.c_str(), maxclustnum, 
                             nmi_tvalue, beta_tvalue, thisfunc);
}
