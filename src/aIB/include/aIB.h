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


#ifndef AIB_H 
#define AIB_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <iostream>
#include <vector>
#include <map>
#include "tools.h"

using namespace std;
using namespace libtools;



struct bestmerge {
  int bm_l;
  int bm_r;
};



class Prm {
 public:
  
  Prm();
  
  Prm(unsigned int xx,unsigned int yy,double Bbeta, int isUniformpriors);
  
  Prm(unsigned int xx,unsigned int yy,double Bbeta, int isUniformpriors, vector <int>);

  Prm(vector < vector <double> > *MM,double Bbeta, int isUniformpriors, vector <int> );

  ~Prm();

  unsigned int getX(); 
  unsigned int getY();
  double getBeta();
  int getUniformPriors();
  const vector <int> * getTsizes();

 protected:
  unsigned int X,Y;
  double beta;
  vector <int>* Tsizes;
  int Uniformpriors;
};


class Inp {
 public:
  Inp(vector < vector <double> > * M, int Uniform);
  ~Inp();
  double getI();
  double getHx();
  double getHy();
  const vector <vector <double> > *get_Py_x();
  const vector <double> *get_Px();

 protected:
  vector < vector <double> > *Py_x;
  vector < vector <double> > *Pxy;
  vector <double> *Py;
  vector <double> *Px;
  MIvalues initialmi;
};


class TmpT{
 public:
  TmpT(Inp * input, Prm * currprm);
  ~TmpT();
  int get_tmpT_size();
  vector<double>* get_Py_t(int ii);
  vector<int>* get_Pt_x();
  double get_Info_Ity(int ii);
  double get_Pt(int ii);
  int get_tmpt_size();
  double get_Info_Ht(int ii);
  void updatemergelog(unsigned int thismerge, int l,int r);
  void updateDelta(unsigned int thismerge,double deltavalue);
  void updateclustering(int thisr, int newval);
  void updatePy_t(int thisl, int thisr);
  void updatePt(int thisl, int thisr);
  void updateInfoValues(int thisl, int thisr,double beta, double thisdelta);
  void updateDeltaL(vector< vector<double> >* delta, TmpT* loctmpt, int thisl, int thisr, double beta,functionals ff);

 protected:
  vector <vector <int> > *MergeLog;
  vector <double>  *Info_Ity;
  vector <double>  *Info_Ht;
  vector <double>  *Info_deltaL;
  //AAACTHUNG tmpt_size is the index of the last element in the vectors/matrix not its size - this mess up is because in the original MATLAB notation vectors start from 1 - so let's stick to Slonim ...
  int tmpt_size;
  int prevTsize;
  double pi_l,pi_r;
  vector <int>  *Pt_x;
  vector < vector <double> >  *Py_t;
  vector <double>  *Pt;
};



class TT_elem {
 public:
  TT_elem(TmpT *r_tmpt, Inp *r_inp, Prm * r_prm);
  double get_Ity();
  vector <int> * get_Pt_x();
  ~TT_elem();

 protected:
  vector <int>  *Pt_x;
  double Ity;
  double Ht;
  double Ity_div_Ixy;
  double Ht_div_Hx;
  int size;
};


#endif
