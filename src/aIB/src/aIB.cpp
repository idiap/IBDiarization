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


#include "aglobal.h"
#include "aIB.h"
#include "fileIO.h"
#include "tools.h"
#include "various_exceptions.h"
#include "aibfunctions.h"

using namespace std;

//various aIB clustering parameters

template <class T> T GetValue (T a) { return a; }

Prm::Prm(): X(0),Y(0),beta(0), Tsizes(new vector<int>(1,1)) , Uniformpriors(0) {}

Prm::Prm(unsigned int xx,unsigned int yy,double Bbeta, int isUniformpriors): 
            X(xx),
            Y(yy),
            beta(1/Bbeta),
            Uniformpriors(isUniformpriors)  {}

Prm::Prm(unsigned int xx,unsigned int yy,
         double Bbeta, int isUniformpriors, 
         vector <int> vv): 
            X(xx),
            Y(yy),
            beta(1/Bbeta),
            Uniformpriors(isUniformpriors)  
{ 
    Tsizes = new vector<int> (vv);  
}

Prm::Prm(vector < vector <double> > *MM,
         double Bbeta, 
         int isUniformpriors, 
         vector <int> vv ): 
            X(MM->size()),
            Y((MM[0]).size()),
            beta(1/Bbeta),
            Uniformpriors(isUniformpriors) 
{ 
    Tsizes = new vector<int> (vv); 
}

Prm::~Prm(){ 
    Tclean_vec(Tsizes);
    delete Tsizes;
}


unsigned int Prm::getX() { return GetValue(X);  }
unsigned int Prm::getY() { return GetValue(Y);  }
double Prm::getBeta() {return GetValue(beta);}
int Prm::getUniformPriors() { return GetValue(Uniformpriors);  }
const vector <int> *  Prm::getTsizes() { return GetValue(Tsizes); }


//various aIB input parameter processing
//Inp::Inp is the only initializer that use and checks that everything is actually in ... 
//     so you can't initialize without the actual values

Inp::Inp(vector < vector <double> > * M, int Uniform)  
{
  
    if (!libtools::check_if_matrix(M)) { 
      throw notamatrix(); 
    }

    int m=M->size();
    int n=(*M)[0].size();

    Py_x = new vector < vector <double> > ( m, vector <double> (n,0) );
    Pxy = new vector < vector <double> > ( m, vector <double> (n,0) );
    Py = new vector <double> (n,0);
    Px = new vector <double> (m,0);

    if ( Uniform ) {
      vector < vector <double> > tmp_matrix (m, vector <double> (n,0) );
      libtools::NormMat(M,&tmp_matrix); 
      libtools::NormMat(M,Py_x);
       
      double totsum=libtools::summatrix(&tmp_matrix); 
      libtools::dividematrixscalar(&tmp_matrix,Pxy,totsum);
      
      tmp_matrix.clear();
    }
    else {
      libtools::NormMat(M,Py_x);
      double totsum=libtools::summatrix(M);
      libtools::dividematrixscalar(M,Pxy,totsum);
    }

    libtools::sumcolumn(Pxy,Px); 
    libtools::sumrow(Pxy,Py);  
    libtools::computeMIvalues(Pxy,&initialmi);

    cout << " The initial MI values are: " << initialmi.I << " " << initialmi.Hx << " " << initialmi.Hy << " \n";
}


Inp::~Inp() { 
    Tclean_vec(Px); delete Px;
    Tclean_vec(Py); delete Py;
    Tclean_mat(Py_x); delete Py_x;
    Tclean_mat(Pxy);  delete Pxy;
}


double Inp::getI() { return GetValue(initialmi.I);  }
double Inp::getHx() { return GetValue(initialmi.Hx);  }
double Inp::getHy() { return GetValue(initialmi.Hy);  }
const vector < vector <double> > * Inp::get_Py_x() { return Py_x;  }
const vector <double> * Inp::get_Px() {return Px;}







//various TmpT input parameter processing
//TmpT::TmpT is the only initializer that use and checks that everything is 
//actually in ... so you can't initialize without the actual values

TmpT::TmpT(Inp* input, Prm* currprm) : tmpt_size(currprm->getX() - 1)
{
    MergeLog = new vector <vector <int> > (2,vector <int> (currprm->getX(),0));
    Info_Ity = new vector <double> (currprm->getX(),0);
    Info_Ht = new vector <double> (currprm->getX(),0);
    Info_deltaL = new vector <double> (currprm->getX(),0);


    (*Info_Ity)[currprm->getX()-1]=input->getI();
    (*Info_Ht)[currprm->getX()-1]=input->getHx();


    const vector < vector <double> > *tmp_point = input->get_Py_x();
    Py_t = new vector < vector <double> > (tmp_point->size(),
                                           vector <double> ((*tmp_point)[0].size(),
                                                             0));
    copymatrix(const_cast<vector < vector <double> > *>(tmp_point),Py_t);

    const vector <double> *tmp_vec=input->get_Px();
    Pt = new vector <double> (tmp_vec->size(),0);
    (*Pt).swap(*(const_cast< vector <double> *>(tmp_vec)));



    Pt_x = new vector <int> (currprm->getX(),0);
    unsigned int ii;
    for(ii=0; ii < currprm->getX(); ii++) { (*Pt_x)[ii] = ii; } 

    prevTsize = pi_l =pi_r = 0;
}


TmpT::~TmpT(){ 
    Tclean_mat(MergeLog); delete MergeLog;
    Tclean_vec(Info_Ity); delete Info_Ity;
    Tclean_vec(Info_Ht);  delete Info_Ht;
    Tclean_vec(Info_deltaL); delete Info_deltaL;
    Tclean_vec(Pt_x); delete Pt_x;
    Tclean_mat(Py_t); delete Py_t;
    Tclean_vec(Pt); delete Pt;
}


int TmpT::get_tmpT_size() { return tmpt_size;  }
vector <double> *  TmpT::get_Py_t(int ii) { return &((*Py_t)[ii]); }
double TmpT::get_Pt(int ii) { return (*Pt)[ii]; }
vector <int> *  TmpT::get_Pt_x() { return Pt_x; }
int TmpT::get_tmpt_size() { return tmpt_size; }
double TmpT::get_Info_Ity(int ii) { return (*Info_Ity)[ii]; }
double TmpT::get_Info_Ht(int ii) { return (*Info_Ht)[ii];  }

void TmpT::updatemergelog(unsigned int thismerge, int l,int r)
{
  if ( thismerge > (*MergeLog)[0].size() ) { throw errorduringmerge(); }
  (*MergeLog)[0][thismerge] = l;
  (*MergeLog)[1][thismerge] = r;

  
  return;
}

void TmpT::updateDelta(unsigned int thismerge,double deltavalue)
{
    if (thismerge > (*Info_deltaL).size()) { 
        throw errorduringmerge();  
    }
    (*Info_deltaL)[thismerge] = deltavalue;
    return;
}


void TmpT::updateclustering(int thisr, int newval)
{ 
    std::replace((*Pt_x).begin(), (*Pt_x).end(), thisr, newval);  
}

void TmpT::updatePy_t(int thisl, int thisr)
{
    int newval=thisl;
    pi_l=(*Pt)[thisl]/((*Pt)[thisl]+(*Pt)[thisr]);
    pi_r=(*Pt)[thisr]/((*Pt)[thisl]+(*Pt)[thisr]);

    //temporary vectors for averaging
    vector <double> tmpvec ((*Py_t)[thisl]);
    unsigned int ii=0;
    for (ii=0; ii < tmpvec.size(); ii++)
    { tmpvec[ii] = pi_l * ((*Py_t)[thisl][ii]) + pi_r * ((*Py_t)[thisr][ii]); }

    std::copy(tmpvec.begin(),tmpvec.end(),(*Py_t)[newval].begin());
    tmpvec.clear();

    vector <double> tmpvec_zeros ((*Py_t)[thisl].size(),0);
    std::copy(tmpvec_zeros.begin(),tmpvec_zeros.end(),(*Py_t)[thisr].begin());
    tmpvec_zeros.clear();

    return;
}


void TmpT::updatePt(int thisl, int thisr)
{
    int newval=thisl;
    (*Pt)[newval]=(*Pt)[thisl]+(*Pt)[thisr];
    (*Pt)[thisr]=0;

    prevTsize = tmpt_size;
    tmpt_size = tmpt_size -1;

    return;
}

void TmpT::updateInfoValues(int thisl, int thisr,double beta, double thisdelta)
{
    double prevL= (*Info_Ity)[prevTsize] - beta*(*Info_Ht)[prevTsize];
    vector <double> tmpvec (2,0); tmpvec[0]=pi_l; tmpvec[1]=pi_r;
    double H = entropy(&tmpvec);
    tmpvec.clear();

    (*Info_Ht)[tmpt_size] = (*Info_Ht)[prevTsize] - (*Pt)[thisl]*H;
    (*Info_Ity)[tmpt_size] = prevL - thisdelta + beta*(*Info_Ht)[tmpt_size];

    return;
}



void TmpT::updateDeltaL(vector< vector<double> > *thisdelta, TmpT * thisloctmpt, int thisl, int thisr, double beta,functionals ff)
{
    unsigned int i;
    for (i=0; i< (*thisdelta).size(); i++) { 
        (*thisdelta)[i][thisr] = numeric_limits<double>::max( );   
    }
    for (i=0; i< (*thisdelta)[0].size(); i++) { 
        (*thisdelta)[thisr][i] = numeric_limits<double>::max( ); 
    }

    omp_set_num_threads(ff.threads_num);
    #pragma omp parallel for
    for (i=0; i< (*thisdelta).size(); i++) {
        if ( (*thisdelta)[i][thisl] <  numeric_limits<double>::max() ) {
            (*thisdelta)[i][thisl] = CalcDelta(thisloctmpt->get_Py_t(i),
                                               thisloctmpt->get_Pt(i), 
                                               thisloctmpt->get_Py_t(thisl),
                                               thisloctmpt->get_Pt(thisl),beta);
        }
    }
    #pragma omp taskwait

    omp_set_num_threads(ff.threads_num);
    #pragma omp parallel for
    for (i=0; i< (*thisdelta)[0].size(); i++) {
        if ((*thisdelta)[thisl][i] < numeric_limits<double>::max( )) {
            (*thisdelta)[thisl][i] = CalcDelta(thisloctmpt->get_Py_t(thisl),
                                               thisloctmpt->get_Pt(thisl), 
                                               thisloctmpt->get_Py_t(i),
                                               thisloctmpt->get_Pt(i),beta);
        }
    }
    #pragma omp taskwait
}


TT_elem::TT_elem(TmpT *r_tmpt, Inp *r_inp, Prm * r_prm)
{
    size = r_tmpt->get_tmpt_size();
    Ity = r_tmpt->get_Info_Ity(r_tmpt->get_tmpt_size());
    Ht = r_tmpt->get_Info_Ht(r_tmpt->get_tmpt_size());
    Ity_div_Ixy = Ity / r_inp->getI();
    Ht_div_Hx = Ht / r_inp->getHx();

    cout << " Size: " << size << " Ity: " << Ity << " Ht: " << Ht 
         << " Ity_div_Itx:  " << Ity_div_Ixy << "  Ht_div_Hx:   " 
         << Ht_div_Hx << "\n";

    vector <int> *clustnames;
    vector<int>::iterator it;
    clustnames = new vector <int> ( *(r_tmpt->get_Pt_x()) ); 
    Pt_x = new vector <int> ( *(r_tmpt->get_Pt_x()) ); 


    sort ((*clustnames).begin(),(*clustnames).end());
    it = unique ((*clustnames).begin(), (*clustnames).end());
    (*clustnames).resize( it - (*clustnames).begin() );

    vector<int>::iterator itlong;
    for (itlong = (*Pt_x).begin(); itlong != (*Pt_x).end(); itlong++)
    {

        *itlong=find((*clustnames).begin(), (*clustnames).end(), *itlong)-(*clustnames).begin();
    }

    delete clustnames;
}

vector <int> * TT_elem::get_Pt_x() { return Pt_x; }
double TT_elem::get_Ity() { return Ity; }


TT_elem::~TT_elem()
{
    Tclean_vec(Pt_x); delete Pt_x;
}

