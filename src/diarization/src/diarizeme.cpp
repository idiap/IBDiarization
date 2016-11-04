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
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>

using namespace std;
namespace po = boost::program_options;

#define PRECISION 1.0e7

#include "global.h"
#include "hmmreadwrite.h"
#include "scpread.h"
#include "featconfig.h"
#include "extract_features.h"

#include "aglobal.h"
#include "aIB.h"
#include "fileIO.h"
#include "tools.h"
#include "aibfunctions.h"

#include "realign.h"

int main (int argc, char** argv)
{
  time_t t1,t2;
  time(&t1);
  double start_time, aibfeat_end_time, aib_end_time, prg_end_time;
  
  cout << " AIB diarization started at " << t1/3600 << " \n";



  std::string ConfigFile, scpFile;
  
  //variables for Deepus input
  std::vector< std::string> mfccFile, tdoaFile, otherFile;
  vector<FeatStream>::iterator featIter;
  
  ConfigVars configVars;
  FeatStream currFeat;
  FeatType currFeatType ;
  string strarg;
  vector<FeatStream> &m_feat = configVars.m_feats;


  try {
    po::options_description generic("aIB Clustering - fabio.valente@idiap.ch ");
    generic.add_options()
      ("version,v","versions 0.00")
      ("help,h","no help for the time being")
      ("config,c",po::value<std::string>(&ConfigFile),"Name of a file of a configuration.")
      ;
    
    po::options_description config("Configuration");
    config.add_options()
      ("recid,d",po::value<std::string>()->default_value("ID")->composing(), "the recording id")
      ("outdir,o",po::value<std::string>()->default_value("./")->composing(),"output directory")
      ("tmpdir,t",po::value<std::string>()->default_value("./")->composing(),"temporary directory")
      ("scp,s",po::value<std::string>(&scpFile)->composing(),"speech/non-speech file")
      ("mfcc",po::value<std::vector<std::string> >(&mfccFile)->multitoken()->composing(),"mfcc feature file and mfcc weight")
      ("tdoa",po::value<std::vector<std::string> >(&tdoaFile)->multitoken()->composing(),"tdoa feture file and tdoa weight")
      ("other",po::value<std::vector<std::string> >(&otherFile)->multitoken()->composing(),"other feature files and other weights")
      ("maxdur",po::value<int>()->default_value(250)->composing(),"Maximum duration of input segments")
      ("maxclust",po::value<int>()->default_value(10)->composing(),"maximum number of allowed clusters")
      ("minHMMdur",po::value<int>()->default_value(250)->composing(),"minimum state duration in the HMM/RKL decoder")
      ("NMIvalue",po::value<double>()->default_value(0.5)->composing(),"normalized mutual information threshold")
      ("beta",po::value<double>()->default_value(10)->composing(),"beta regularizer")
      ("nthread",po::value<int>()->default_value(4)->composing(),"number of threads")
      ;
    
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);
    
    po::options_description config_file_options;
    config_file_options.add(config);
    
    po::options_description visible("Options");
    visible.add(generic).add(config);
    
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    
    notify(vm);
    
    if (vm.count("help")) {
      cout << visible << "\n";
      return 0;
    }
    
    if (vm.count("version")) {
      cout << "aIB feature extraction - C++ code - version 0.00\n";
      return 0;
    }
    
    
    if (vm.count("config")) {
        ifstream ifs(ConfigFile.c_str());
        if(ifs) { 
            store(parse_config_file(ifs, config_file_options), vm); 
            notify(vm); 
        }
        else {
          cout << " Cannot open the config file " <<  ConfigFile.c_str() << "\n";   return 0;
        }
    }
    
    if (!vm.count("scp")) {
        cout << "Scp file not found!!: Exiting \n"; 
        return 0;
    }
    else {
        ifstream ifs(scpFile.c_str());
        if (!ifs) {  
            cout << " Cannot open the scp file " <<  scpFile.c_str() << "\n"; 
            return 0; 
        }
        else { 
            configVars.scpfile = scpFile.c_str();  
        }
    }
    
    cout << "\n\n";

    if (vm.count("mfcc")) {
        currFeatType = MFCC;
        if (mfccFile.size() < 2) { 
            cout << "invalid passing to mfcc "; 
            return 0; 
        }
        currFeat.m_feattype=currFeatType;
        currFeat.m_filename = mfccFile[0].c_str();
        currFeat.m_wt = atof(mfccFile[1].c_str());
        cout << " MFCC:   " << currFeat.m_filename << "  weight MFCC  " << currFeat.m_wt << " \n";
        m_feat.push_back(currFeat);
    }

    if (vm.count("tdoa")) {
        currFeatType = TDOA;
        if (tdoaFile.size() < 2) { 
            cout << "invalid passing to tdoa features "; 
            return 0; 
        }
        currFeat.m_feattype=currFeatType;
        currFeat.m_filename = tdoaFile[0].c_str();
        currFeat.m_wt = atof(tdoaFile[1].c_str());
        cout << " TDOA:   " << currFeat.m_filename << "  weight TDOA " << currFeat.m_wt << " \n";
        m_feat.push_back(currFeat);
    }


    if (vm.count("other")) {
        currFeatType = OTHER;
        vector<string> strs;
        if (otherFile.size() < 2) {
            cout << "invalid passing to other features"; 
            return 0; 
        }
        
        for (unsigned int ii=0; ii < otherFile.size(); ii = ii + 2)
          {
            currFeat.m_feattype=currFeatType;
            currFeat.m_filename = otherFile[ii].c_str();
            currFeat.m_wt = atof( otherFile[ii+1].c_str() );
            cout << " OTHER:   " << currFeat.m_filename << "  weight OTHER " << currFeat.m_wt << " \n";
            m_feat.push_back(currFeat);
          }
    }

    
    if (m_feat.size() == 0 ){
          cerr << "No feature value set\n"
               << "please set with --mfcc/--tdoa/--other <filename> <feat_wt>"
               << endl;
          exit(1);
    }

    vector<FeatStream>::iterator itfea;
    float sum_w=0;
    for(itfea=m_feat.begin(); itfea < m_feat.end(); itfea++) {
        if (itfea->m_wt < 0) {
            cout << "Negative feature weights !!"; 
            return 2;
        }
        sum_w+= itfea->m_wt;
    }
    
    if (sum_w != 1) {
        cout << "Sum of feature weights is not one !!"; 
        return 0; 
    }

    if (vm.count("maxdur")) {
       configVars.m_maxDur=vm["maxdur"].as<int>();
    }

    if (vm.count("minHMMdur")) { 
        configVars.m_minDur=vm["minHMMdur"].as<int>();
    }

    cout << "\n Maximum Segment Duration " << configVars.m_maxDur << " \n";

    configVars.m_doViterbi=1;

    
    if (vm.count("recid")) {
        configVars.id = vm["recid"].as<string>();  
    }
    
    if (vm.count("tmpdir")) { 
        configVars.m_tmpDir = vm["tmpdir"].as<string>(); 
    }

    string outfile_s(vm["outdir"].as<string>());
    outfile_s.append("/");
    outfile_s.append(vm["recid"].as<string>());
    outfile_s.append(".clust.out");
    const char *outfile=outfile_s.c_str();
    
    int maxclustnum = (int)vm["maxclust"].as<int>();
    configVars.m_maxClusters = maxclustnum;
    double nmi_tvalue = (double)vm["NMIvalue"].as<double>();
    double beta_tvalue = (double)vm["beta"].as<double>();
    
    functionals thisfunc;
    thisfunc.threads_num = (int)vm["nthread"].as<int>();
    
    if ( maxclustnum < 1 ) {
        cout << "The maximum number of cluster must be larger then 1\n";
        return 0;
    }
    
    if (nmi_tvalue <=0 || nmi_tvalue >1) {
        cout << "The NMI threshold must be included between 0 and 1\n";
        return 0;
    }
    
    if (beta_tvalue <= 0) {
        cout << "The beta value must be striclty larger then zero \n";
        return 0;
    }
    
    if (thisfunc.threads_num <= 0) {
        cout << "The number of threads must be larger then zero";
        return 0;
    }
    
    cout << "\n\n AIB is running with the following parameters\n\n" ;
    cout << " Maximum number of clusters possible: \t\t" << maxclustnum << "\n";
    cout << " Normalized Mutual Information threshold: \t\t" << nmi_tvalue << "\n";
    
    cout << " Beta value: \t\t" << beta_tvalue << "\n";
    cout << " Number of threads: \t\t" << thisfunc.threads_num << "\n";
    
    ///This is the diarization core
    vector <vector <float> > postmat;

    cout << "\n\n";

    start_time = get_cpu_time();
    compute_features(configVars,postmat);
    aibfeat_end_time = get_cpu_time();
    cout << "compute features after: "<< aibfeat_end_time - start_time << endl;
    cout << "\n\n";
    
    load_prepare_and_cluster(postmat,outfile,maxclustnum,nmi_tvalue,beta_tvalue,thisfunc);
    aib_end_time = get_cpu_time();
    
    cout << "\n\n";
    
    dorealign(configVars);
    prg_end_time = get_cpu_time();
    
    time(&t2);
    cout << "Diarization stopped at " << t2/3600 << "  \n";
    cout << "Running time " << difftime (t2,t1) << " seconds  \n";
    cout << "compute features after: "<< aibfeat_end_time - start_time << endl;
    cout << "aib clustering after: " << aib_end_time - aibfeat_end_time << endl;
    cout << "realignment after: " << prg_end_time - aib_end_time << endl;
    
  } catch(std::exception& e) {
     cout << "Something went wrong : " << e.what() << "\n";
     return 1;
  }
}
