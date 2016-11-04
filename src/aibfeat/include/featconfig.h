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


#ifndef FEATCONFIG_H
#define FEATCONFIG_H


enum FeatType {MFCC, PLP, TDOA, OTHER};
class FeatStream {
   public:
      FeatType m_feattype;
      string m_filename ;
      float m_wt;
      FILE *featinp; 
};

class ConfigVars{
   public:
      vector<FeatStream> m_feats ; 
      int m_maxClusters;
      int m_maxDur; 
      int m_modelSelParam;
      bool m_doViterbi;
      int m_minDur;
      string m_tmpDir;
      string scpfile;
      string id;
};


#endif
