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


#ifndef LIBTOOLS
#define LIBTOOLS
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <functional>
#include <numeric>
#include <iostream>
#include <iterator>
#include <sstream>
#include <fstream>

//#include "/idiap/home/shyella/boost_1_45_0/boost/smart_ptr.hpp"
#include "boost/smart_ptr.hpp"
#define TOL 1E-7

//using namespace boost;
using namespace std;
using namespace boost;



template <class T> void Tprint (vector < vector<T> >*a)
    {
      typename  vector <vector< T > >::iterator thisIterator;
      ostream_iterator< T > out_it (cout," ");
      for (thisIterator=(*a).begin(); thisIterator!=(*a).end(); thisIterator++)
	{ copy((*thisIterator).begin(),(*thisIterator).end(),out_it);  cout << "\n"; }
      return;
    }


template <class T> void Tprint (vector<T> *a)
    {
      ostream_iterator< T > out_it (std::cout," ");
      copy((*a).begin(),(*a).end(),out_it);
      cout << "\n";
      return;
    }



template <class T> void Tprint (vector < vector<T> >*a, const char* outfilename)
{
  typename  vector <vector< T > >::iterator thisIterator;
  ofstream myfile (outfilename);
  ostream_iterator< T > out_it (myfile," ");
  for (thisIterator=(*a).begin(); thisIterator!=(*a).end(); thisIterator++)
    { copy((*thisIterator).begin(),(*thisIterator).end(),out_it);  myfile << "\n"; }
  myfile.close();
  return;
}

template <class T> void Tprint (vector<T> *a, const char* outfilename)
    {
      ofstream myfile (outfilename);
      ostream_iterator< T > out_it (myfile," ");
      copy((*a).begin(),(*a).end(),out_it);
      myfile << "\n";
      myfile.close();
      return;
    }



template <class T> void Tclean_vec (vector<T> *a)
{
  (*a).clear();
  vector<T>().swap(*a);
}

template <class T> void Tclean_mat (vector <vector< T > >*a)
{
  typename  vector <vector< T > >::iterator deleteIterator;
  for (deleteIterator= (*a).begin(); deleteIterator!= (*a).end();  deleteIterator++)
    {
      Tclean_vec(&(*deleteIterator));
    }

     
    (*a).clear();
    vector <vector<T> >().swap(*a);
}



namespace libtools
{

  struct MIvalues {
    double I;
    double Hx;
    double Hy;
  };


  double MyLog(double *) ;
  inline double MyLog(double  ) ;
  inline int checkprobvector(vector < double > *vv );
  int checkprobvector(boost::shared_ptr< vector<double> > *vv );
  inline double kldivergence(vector <double> *v1, vector <double> *v2, int Chk_q);
  double jsdivergence(vector <double> *v1, vector <double> *v2, double pi_1, double pi_2);
  double entropy(vector <double> *v1);
  int check_if_matrix(vector< vector< double> > *mm);
  int check_jointmatrix(vector< vector< double> > *mm);
  void computeMIvalues(vector< vector< double> > *mm, MIvalues *val);
  void NormMat(vector< vector< double> > *mm_in, vector< vector< double> > *mm_out);
  double summatrix(vector< vector< double> > *mm);
  void printmatrix(vector< vector< double> > *mm);
  void printvector(vector< double> *vv);
  void dividematrixscalar(vector< vector< double> > *mm_in, vector< vector< double> > *mm_out,double scalar);
  void sumcolumn(vector< vector< double> > *mm_in, vector< double>  *vec_out);
  void sumrow(vector< vector< double> > *mm_in, vector< double>  *vec_out);
  void copymatrix(vector< vector< double> > *mm_in, vector< vector< double> > *mm_out);
  void remove_emptycolumn_fromcounts(vector< vector<double> > *mm_in,vector< vector<double> > *mm_out);
  double get_cpu_time();
  

  

}

#endif
