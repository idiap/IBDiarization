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


#ifndef AIBFUNCTIONS_H 
#define AIBFUNCTIONS_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <iostream>
#include <vector>
#include <map>
#include "tools.h"
#include "aglobal.h"

using namespace std;
using namespace libtools;




void InitDeltaL(vector< vector<double> > *,TmpT* actualTmpT,Prm* actualPrm, functionals ff);
double CalcDelta(vector<double>* P_ytl, double p_tl, vector<double>* P_ytr, double p_tr, double beta);
void Find_Bestmerge(vector< vector<double> > *, bestmerge *);
void MergeAndUpdate(vector< vector<double> > *delta, TmpT * loctmpt, int thisl, int thisr, int thismerge, Prm* thisprm, functionals ff);
int check_if_sholdsave(Prm *lprm,int thisiter);

typedef boost::shared_ptr<TT_elem> TT_elem_ptr;
void DoAIBclustering(vector <vector <double> >* InputM,double beta,int Uprior,vector <int> listofrecords,double t_value_nmi,vector <int>* clustering_out_vec, functionals ff);
void compute_modelselection_values(map<int,TT_elem_ptr>*history,vector <int> *clustersol, double Ixy, double nmithreshold);
void coolprint(int a,int b);
void  store_misolutions(map<int,TT_elem_ptr>* history, vector <double> *sol, double Ixy);

void load_prepare_and_cluster(const char* inputmatrix, const char* outputfile, unsigned int maxclustnum, double nmi_tvalue, double beta_tvalue,functionals ff);
void load_prepare_and_cluster(vector <vector <float> >&II , const char* outputfile,unsigned int maxclustnum, double nmi_tvalue, double beta_tvalue, functionals ff);
#endif
