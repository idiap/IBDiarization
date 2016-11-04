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
#include "sys/time.h"

#include "tools.h"
#include "various_exceptions.h"
#include "boost/config.hpp"
#include "boost/smart_ptr.hpp"

using namespace boost;
using namespace std;


namespace libtools
{

  double  MyLog(double* value) 
  {

    if (*value <= 0) {throw log_negative();}
    else { return log(*value); }
    
  }

  double  MyLog(double value) 
  {

    if (value <= 0) {throw log_negative();}
    else { return log(value); }
    
  }


  int checkprobvector (vector<double> *vv )
  {
    unsigned int isprobvector = 1;  
    vector <double>::iterator iter;
    for (iter = vv->begin(); iter != vv->end(); ++iter)
      { if ( *iter < 0) { isprobvector = 0;  }  }

    
    double vsum = 0;
    vsum=std::accumulate( vv->begin(), vv->end(), 0.0E-20 );
    if ( fabs( vsum - (double)1.0 ) > TOL ) { isprobvector = 0;  }


    //cout << "sum " << fabs( vsum - (double)1.0 ) << " " <<  (*vv).size()<<"\n";
    return isprobvector;

  }



  double kldivergence(vector <double> *v1, vector <double> *v2, int Chk_q)
  {
    if ( !(checkprobvector(v1)) ) { throw notaprobvector(); }
    if ( !(checkprobvector(v2)) ) { throw notaprobvector(); }
    
    double klvalue=0.0E-20;
    
    vector<double>::iterator posmin;
    posmin = min_element (v2->begin(), v2->end());
    
    if (v1->size() != v2->size() ) {cout << "not a valid divergence\n"; exit(1);}
    
    vector<double>::iterator kliterator_v1;
    vector<double>::iterator kliterator_v2;
    
    
    for (kliterator_v1 = v1->begin(), 
	   kliterator_v2 = v2->begin(); kliterator_v1 != v1->end(), 
	   kliterator_v2 != v2->end(); ++kliterator_v1,++kliterator_v2 )
      {
	
	if (*kliterator_v1 > 0)
	  //if (*kliterator_v1 > 0 && *kliterator_v2 > 0)
	  {
	    double tmp=(*kliterator_v1)/(*kliterator_v2);
	    klvalue = klvalue + (*kliterator_v1)* MyLog(tmp);
	    //cout << "KL value " << klvalue << " " << *kliterator_v1 << " " << *kliterator_v2 << " " << MyLog(tmp) <<"\n";
	  }
      }
    
    if ( !(klvalue >=0) ){ cout << "problem with a KL " << klvalue  << "\n\n"; exit(1);}
    
    return klvalue;
  }

  

  double jsdivergence(vector <double> *v1, vector <double> *v2, double pi_1, double pi_2)
  {
    
    double jsvalue=0;

    if ( !(checkprobvector(v1)) ) { throw notaprobvector(); }
    if ( !(checkprobvector(v2)) ) { throw notaprobvector(); }

    if ( fabs(pi_1+pi_2 - 1) > TOL || pi_1 < 0 || pi_2 <0 ) 
      { throw notaprobvector(); }

    vector <double> average ((*v1).size(),0);
    unsigned int i;
    for(i=0; i< (*v1).size(); i++)
      { average[i]=( (*v1)[i]*pi_1 + (*v2)[i]*pi_2 );  }

    
    
    double kl_1;
    if (pi_1 > 0) {  kl_1=kldivergence(v1,&average,0);} else {kl_1=0;}
    
    double kl_2;
    if (pi_2 > 0) { kl_2=kldivergence(v2,&average,0);} else {kl_2=0;}
    
    jsvalue= pi_1*kl_1 + pi_2*kl_2;

    return jsvalue;
  }





  double function_plogp(double a)
  {
    if (a ==0) return 0;
    else {return -a*MyLog(a); }
  }



  double entropy(vector <double> *v1)
  {
    if ( !(checkprobvector(v1)) ) { throw notaprobvector(); }

    vector <double>::iterator v1_iterator;
    vector <double> plogp ((*v1).size());
    
    transform ((*v1).begin(), (*v1).end(), plogp.begin(),function_plogp);
    double totentropy = std::accumulate( plogp.begin(), plogp.end(), 0.0E-20 );
    return totentropy;
  }


  int check_if_matrix(vector< vector< double> > *mm)
  {
    int ismatrix=1;
    vector < vector<double> >::iterator rowiter;
    
    if (mm->size() < 2) {  ismatrix=0; }

    unsigned int columndim = (unsigned int) (*mm)[0].size();
        
    for (rowiter = mm->begin(); rowiter != mm->end(); ++rowiter)
      { if ((*rowiter).size() != columndim) { ismatrix=0; }  }

    return ismatrix;
  }


  int check_jointmatrix(vector< vector< double> > *mm)
  {
    double sumpxy=0;
    int isjointmatrix=1;

    if (!check_if_matrix(mm)) {throw notamatrix();}
    
    vector < vector<double> >::iterator rowiter;
    for (rowiter = mm->begin(); rowiter != mm->end(); ++rowiter)
      { 
	vector<double> *tt = &(*rowiter);
	sumpxy += std::accumulate( tt->begin(), tt->end(), 0.0E-20 );
      }
    

    if ( fabs(sumpxy-1) > TOL) { isjointmatrix=0; }

    vector <double>::iterator columniter;
    for (rowiter = mm->begin(); rowiter != mm->end(); ++rowiter)
      for (columniter=(*rowiter).begin(); columniter != (*rowiter).end(); ++columniter)
	{  if (*columniter < 0) { isjointmatrix=0; }  }

    return isjointmatrix;
  }




  void NormMat(vector< vector< double> > *mm_in, vector< vector< double> > *mm_out)
  {
     if (!check_if_matrix(mm_in)) {throw notamatrix();}
     if (!check_if_matrix(mm_out)) {throw notamatrix();}

     if (mm_in->size() != mm_out->size()) {throw dimismatch();}
     if ((*mm_in)[0].size() != (*mm_out)[0].size()) {throw dimismatch();}
     
     unsigned int m=mm_in->size();
     unsigned int n=(*mm_in)[0].size();


     double sums[mm_in->size()];
     unsigned int i=0;
     vector < vector<double> >::iterator rowiter_in;
     for (rowiter_in = mm_in->begin(); rowiter_in != mm_in->end(); ++rowiter_in, ++i)
      {
	sums[i] = std::accumulate( (*rowiter_in).begin(), (*rowiter_in).end(), 0.0E-20 );
      }

     unsigned int p,q =0;
     for (p=0; p < m; p++)
       for (q=0; q < n; q++)
	 {
	   if (sums[p] > 0) { (*mm_out)[p][q] = ((*mm_in)[p][q]) / sums[p]; }
	   else { (*mm_out)[p][q]=0; }
	 }
  }


  void dividematrixscalar(vector< vector< double> > *mm_in, vector< vector< double> > *mm_out,double scalar)
  {
    if (!check_if_matrix(mm_in)) {throw notamatrix();}
    if (!check_if_matrix(mm_out)) {throw notamatrix();}

    if (mm_in->size() != mm_out->size()) {throw dimismatch();}
    if ((*mm_in)[0].size() != (*mm_out)[0].size()) {throw dimismatch();}
    
    unsigned int m=mm_in->size();
    unsigned int n=(*mm_in)[0].size();
    
    unsigned int p,q;
      for(p=0; p < m; p++)
	for (q=0; q<n; q++)
	  {  (*mm_out)[p][q] = (*mm_in)[p][q] / scalar; }

  }




  void copymatrix(vector< vector< double> > *mm_in, vector< vector< double> > *mm_out)
  {
    if (!check_if_matrix(mm_in)) {throw notamatrix();}
    if (!check_if_matrix(mm_out)) {throw notamatrix();}

    if (mm_in->size() != mm_out->size()) {throw dimismatch();}
    if ((*mm_in)[0].size() != (*mm_out)[0].size()) {throw dimismatch();}
    
    unsigned int m=mm_in->size();
    unsigned int n=(*mm_in)[0].size();
    
    unsigned int p,q;
      for(p=0; p < m; p++)
	for (q=0; q<n; q++)
	  {  (*mm_out)[p][q] = (*mm_in)[p][q]; }
    
  }



  void sumcolumn(vector< vector< double> > *mm_in, vector< double>  *vec_out)
  {
    if (!check_if_matrix(mm_in)) {throw notamatrix();}
    unsigned int m=mm_in->size();
    unsigned int n=(*mm_in)[0].size();
    
    unsigned int p,q;
    double tmpsum=0;
    for(p=0; p < m; p++)
    {   
      tmpsum=0; for (q=0; q<n; q++) { tmpsum+=(*mm_in)[p][q];}
      if (tmpsum==0) 
	{throw emptyrowcolumn();}
      else
	(*vec_out)[p]=tmpsum;
    }  
  }


  void sumrow(vector< vector< double> > *mm_in, vector< double>  *vec_out)
  {
    if (!check_if_matrix(mm_in)) {throw notamatrix();}
    unsigned int m=mm_in->size();
    unsigned int n=(*mm_in)[0].size();
    
    unsigned int p,q;
    double tmpsum=0;
    for (q=0; q<n; q++)
    {  
      tmpsum=0;  
      for(p=0; p < m; p++) 
	{ tmpsum+=(*mm_in)[p][q]; }
      
      if (tmpsum==0)
	{throw emptyrowcolumn();}
      else
	(*vec_out)[q]=tmpsum;
    }

    
  }



  void printmatrix(vector< vector< double> > *mm)
  {
    vector < vector<double> >::iterator rowiter;
    vector <double>::iterator columniter;
    
    for (rowiter = mm->begin(); rowiter != mm->end(); ++rowiter)
      {
	for (columniter=(*rowiter).begin(); columniter != (*rowiter).end(); ++columniter)
	{ cout << "\t" << *columniter;}
	cout << "\n";
      }
  }

  void printvector(vector< double> *vv)
  {
    vector<double>::iterator rowiter;
    
    
    for (rowiter = vv->begin(); rowiter != vv->end(); ++rowiter)
      {	 cout << "\t" << *rowiter;}
    cout << "\n";
  }





  double summatrix(vector< vector< double> >*mm)
  {
    double sum = 0.0E-20;
 
    if (!check_if_matrix(mm)) { throw notamatrix(); }
    vector < vector<double> >::iterator rowiter_in;
     for (rowiter_in = mm->begin(); rowiter_in != mm->end(); ++rowiter_in)
      {
	sum += std::accumulate( (*rowiter_in).begin(), (*rowiter_in).end(), 0.0E-20 );
	//cout << "Sums: " << sum << "\n";
      }
  
     return sum;
  }

  


  void computeMIvalues(vector< vector< double> > *mm, MIvalues *val)
  {
    if (!check_jointmatrix(mm)) {throw notajointmatrix();}
    
    int m=mm->size();
    int n=(*mm)[0].size();

    vector <double> Px (m,0.0E-20);
    int p,q;
    for(p=0; p < m; p++)
      for (q=0; q<n; q++)
	{  Px[p]+= (*mm)[p][q]; }
    
    //for(p=0; p < m; p++) {cout << "Px " << Px[p] << "\n";}
    
    (*val).Hx = entropy(&Px); 
    //cout << "Hx " << (*val).Hx << "\n";

    vector <double> Py (n,0.0E-20);
    for (q=0; q<n; q++)
      for(p=0; p < m; p++)
	{ Py[q] += (*mm)[p][q];  }

    //for(q=0; q < n; q++) {cout << "Py " << Py[q] << "\n";}
    
    (*val).Hy = entropy(&Py); 
    //cout << "Hy " << (*val).Hy << "\n";

    double HXY=0.0E-20;
    for (q=0; q<n; q++)
      for(p=0; p < m; p++)
	{
	  if ((*mm)[p][q] > 0)
	    { HXY = HXY + (*mm)[p][q] * MyLog((*mm)[p][q]);}
	}

    //cout << "TT  " << (*mm)[11][11] * MyLog((*mm)[11][11]) << "  \n";


    //cout << "HXY " << HXY << " \n";
    (*val).I=(*val).Hy + (*val).Hx + HXY;

  }


  void remove_emptycolumn_fromcounts(vector< vector<double> > *mm_in, vector< vector<double> > *mm_out)
  {
    unsigned int m=mm_in->size();
    unsigned int n=(*mm_in)[0].size();
    
    
    vector <double> aa;
    unsigned int p,q;
    double tmpsum=0;
    for (q=0; q<n; q++)
    {   
      tmpsum=0; for (p=0; p < m; p++) { tmpsum+=(*mm_in)[p][q];}
      if (tmpsum==0)  { aa.push_back(q); }
    }  

    
    (*mm_out).resize((*mm_in).size());
    for(vector< vector <double> >::iterator it=(*mm_out).begin();  it!=(*mm_out).end(); it++)
      {  (*it).resize(n - aa.size()); }

    unsigned int idx=0;
    for(q=0; q<n; q++)
      {
	if ( find(aa.begin(),aa.end(),q) == aa.end() ) 
	  {
	    //ACTHUNG this is for avoiding precision problems...
	    for (p=0; p < m; p++)
	      { (*mm_out)[p][idx]=(*mm_in)[p][q]+1E-12; }
	    idx++;
	  }
      }

  }

  double get_cpu_time(){
      return (double)clock() / CLOCKS_PER_SEC;
  }
}

