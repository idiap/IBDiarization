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



int check_if_sholdsave(Prm *lprm,int thisiter)
{
    int foundit=0;

    const vector <int> *list_possible_solutions = lprm->getTsizes();
    vector <int>::const_iterator it;

    for (it=(*list_possible_solutions).begin(); 
           it !=(*list_possible_solutions).end(); it++) { 
        if (*it == thisiter) 
            foundit=1;   
    }
    return foundit;
}



double CalcDelta(vector<double>* P_ytl, double p_tl, vector<double>* P_ytr, double p_tr, double beta)
{
    double Pnew=p_tl+p_tr;
    double pi_l=p_tl/Pnew; double pi_r=p_tr/Pnew; 

    double term1=Pnew * jsdivergence( P_ytl, P_ytr,pi_l,pi_r);

    double term2=0;
    if (beta > 0)
    {
        vector <double> t_prob (2,0);
        t_prob[0]=pi_l; t_prob[1]=pi_r;
        term2=Pnew *entropy(&t_prob);
    }
    return (term1 -beta*term2);
}


void InitDeltaL(vector< vector<double> > *distances,TmpT* actualTmpT,Prm* actualPrm,functionals ff)
{
    cout << " Initializing the Delta -- Computing " 
         << actualTmpT->get_tmpT_size()+1 << " times " 
         << actualTmpT->get_tmpT_size()+1 << " Distances \n";

    unsigned int i=0,j=0;

    time_t t1,t2;
    time(&t1);

    unsigned int ll=actualTmpT->get_tmpT_size();
    unsigned int ll2=actualTmpT->get_tmpT_size();

    omp_set_num_threads(ff.threads_num);

    #pragma omp parallel for firstprivate(j)
    for (i=0; i <= ll; i++)
    {
        for (j=i+1; j<=ll2; j++)
        {
            (*distances)[i][j] = CalcDelta(actualTmpT->get_Py_t(i),
                                           actualTmpT->get_Pt(i),
                                           actualTmpT->get_Py_t(j),
                                           actualTmpT->get_Pt(j),
                                           actualPrm->getBeta());
        }
    }
    #pragma omp taskwait

    time(&t2);
    cout << "Running time " << difftime(t2,t1) << "  \n"; 
    return;
}




void Find_Bestmerge(vector< vector<double> > *dist_matrix, bestmerge *bestmcontainer )
{
    double min=1E30;
    int min_l=-1;
    int min_r=-1;

    unsigned int m = dist_matrix->size();
    unsigned int n= (*dist_matrix)[0].size();

    unsigned int p,q;
    for (p=0; p < m; p++) {
        for (q=0; q < n; q++) {
            if ((*dist_matrix)[p][q] < min) {
                min_l=p; min_r=q; min=(*dist_matrix)[p][q];   
            }
        }
    }

    (*bestmcontainer).bm_l=min_l;
    (*bestmcontainer).bm_r=min_r;
    if (min_l > min_r) { throw somethingwronginclustering(); }
}


void MergeAndUpdate(vector< vector<double> > *delta, TmpT * loctmpt, int thisl, 
                    int thisr, int thismerge, Prm* thisprm, functionals ff)
{
    loctmpt->updatemergelog(thismerge,thisl,thisr); 
    loctmpt->updateDelta(thismerge,(*delta)[thisl][thisr]);

    int newval=thisl;
    loctmpt->updateclustering(thisr,newval);
    loctmpt->updatePy_t(thisl,thisr);
    loctmpt->updatePt(thisl,thisr);
    loctmpt->updateInfoValues(thisl,thisr,thisprm->getBeta(),
                              (*delta)[thisl][thisr]);
    loctmpt->updateDeltaL(delta,loctmpt,thisl,thisr,thisprm->getBeta(),ff);

    return;
}


void DoAIBclustering(vector <vector <double> >* InputM,double beta,int Uprior,
                     vector <int> listofrecords, double t_value_nmi, 
                     vector <int>* clustering_out_vec, functionals ff)
{
  Prm currentprm(InputM,beta,Uprior,listofrecords);
  Inp input(InputM,1);
  TmpT tmptaib(&input,&currentprm);

  vector< vector<double> > DeltaL(tmptaib.get_tmpT_size()+1, 
                                  vector<double>(tmptaib.get_tmpT_size()+1, 
                                                 numeric_limits<double>::max()));
  InitDeltaL(&DeltaL,&tmptaib,&currentprm, ff);
   
  map<int,TT_elem_ptr> TrackSolutions;
  

  unsigned int mergenumcounter=0;
  bestmerge thisbestmerge;
  cout << "\n\n Now Clustering \n\n";
  
  for (mergenumcounter=1; 
          mergenumcounter < currentprm.getX(); 
          mergenumcounter++) {
      coolprint(mergenumcounter,currentprm.getX());
      
      if (check_if_sholdsave(&currentprm , tmptaib.get_tmpT_size()) >0 ) {
          TrackSolutions.insert(std::pair<int,TT_elem_ptr>(tmptaib.get_tmpT_size(), 
                                TT_elem_ptr (new TT_elem(&tmptaib,
                                             &input,
                                             &currentprm))));
	  }

      Find_Bestmerge(&DeltaL,&thisbestmerge);
      MergeAndUpdate(&DeltaL,&tmptaib,thisbestmerge.bm_l,thisbestmerge.bm_r,
                     mergenumcounter-1, &currentprm,ff);
      
  }
  cout << "\n";
  
  //Do the last update
  TrackSolutions.insert(std::pair<int,TT_elem_ptr>(tmptaib.get_tmpT_size(),
                        TT_elem_ptr (new TT_elem(&tmptaib,&input,&currentprm))));
  cout << "Saving this solution \n";
  
  vector <int> clusteringsolution (1,0);
  compute_modelselection_values(&TrackSolutions, &clusteringsolution, input.getI(), t_value_nmi);

  (*clustering_out_vec).resize(clusteringsolution.size());
  copy(clusteringsolution.begin(),clusteringsolution.end(),(*clustering_out_vec).begin());

  /// little clean up
  Tclean_mat(&DeltaL);
}


void  store_misolutions(map<int,TT_elem_ptr>* history, vector <double> *sol, double Ixy)
{
    for( map<int, TT_elem_ptr>::iterator ii=(*history).begin(); ii!=(*history).end(); ++ii) {
        double nmi=(*ii).second->get_Ity() / Ixy;
        (*sol).push_back(nmi);
    }
}


void compute_modelselection_values(map<int,TT_elem_ptr>* history,vector <int> *sol_vec, double Ixy, double nmithreshold)
{
    int found_solution=0;
    int key_to_solution=-1;

    for(map<int, TT_elem_ptr>::iterator ii=(*history).begin(); ii!=(*history).end(); ++ii) {
        double nmi=(*ii).second->get_Ity() / Ixy;
        cout << (*ii).first << "  " << (*ii).second->get_Ity() << "  " << Ixy 
             << "  " << nmi << "  \n";

        if(nmi > nmithreshold && found_solution==0)
        {
            key_to_solution=(*ii).first;
            // found_solution=1;
        }
    }


    //if the nmi threshold is too high save the last possible solution and that's it...
    if(found_solution == 0) { 
        key_to_solution=(*history).end()->first;
    }

    if(key_to_solution ==0) { 
        key_to_solution=1; 
    }

    cout << "Key to solution " << key_to_solution  << " NMI value"  << 
            nmithreshold << "  \n"; 

    map<int, TT_elem_ptr>::iterator itt;
    itt = (*history).find(key_to_solution-1);
    vector <int> *tmppoint = (*itt).second->get_Pt_x();  

    (*sol_vec).resize((*tmppoint).size());
    copy((*tmppoint).begin(),(*tmppoint).end(),(*sol_vec).begin());

    return;
}


void coolprint(int a, int b)
{
    double percentage=1000*((double)a/(double)b);
    if ((int)floor(percentage)%10 == 0) { cerr << "-"; }
}


void load_prepare_and_cluster(const char* inputmatrix, const char* outputfile, 
                              unsigned int maxclustnum, double nmi_tvalue, 
                              double beta_tvalue, functionals ff)
{
    cout << "\n\n Loading the data matrix : "<< inputmatrix << "\n\n";
    vector< vector<double> > AA(1,vector<double>(1,0));
    vector< vector<double> > nozeros_AA(1,vector<double>(1,0));
    readmatrixfromfile_and_resize(inputmatrix,&AA);

    if (libtools::check_if_matrix(&AA)) { 
        cout << " The matrix seems ok\n"; 
    } 
    else  {
        cout << "Error in reading the matrix\n"; 
        throw somethingwrongininputmatrix();
    }

    libtools::remove_emptycolumn_fromcounts(&AA,&nozeros_AA);

    vector<int> v_maxclust(maxclustnum);
    for(unsigned int ii=0; ii<maxclustnum; ii++) { v_maxclust[ii]=ii; }
    vector<int> v_sol(1,0);

    cout << " Final Matrix size after zero removal " << nozeros_AA.size() << " " 
         << nozeros_AA[0].size() << "  \n";

    DoAIBclustering(&nozeros_AA,beta_tvalue,1,v_maxclust,nmi_tvalue,&v_sol,ff);

    cout << " The clustering has finished. Now saving the solution in " 
         << outputfile << endl;

    ///increase the index value to one to match the matlab output

    for (vector<int>::iterator ii=v_sol.begin(); ii !=v_sol.end(); ii++) { 
        *ii=*ii+1;  
    }
    Tprint(&v_sol,outputfile);

    cout << " ...and saying goodbye \n";
}


void load_prepare_and_cluster(vector <vector <float> >&II , 
                              const char* outputfile, unsigned int maxclustnum, 
                              double nmi_tvalue, double beta_tvalue, 
                              functionals ff)
{
    vector< vector<double> > AA(II.size(),vector<double>(II[0].size(),0));

    for (unsigned int i=0; i < II.size(); i++) {  
        copy(II[i].begin(),II[i].end(),AA[i].begin()); 
    }


    vector< vector<double> > nozeros_AA(1,vector<double>(1,0));
    if (libtools::check_if_matrix(&AA)) { 
        cout << " The matrix seems ok\n"; 
    }
    else {
        cout << "Error in reading the matrix\n"; 
        throw somethingwrongininputmatrix();
    }

    libtools::remove_emptycolumn_fromcounts(&AA,&nozeros_AA);

    vector<int> v_maxclust(maxclustnum);
    for(unsigned int ii=0; ii<maxclustnum; ii++) { 
        v_maxclust[ii]=ii; 
    }
    vector<int> v_sol(1,0);

    cout << " Final Matrix size after zero removal " << nozeros_AA.size() 
         << " " << nozeros_AA[0].size() << "  \n";

    DoAIBclustering(&nozeros_AA,beta_tvalue,1,v_maxclust,nmi_tvalue,&v_sol,ff);

    cout << " The clustering has finished. Now saving the solution in " 
         << outputfile << endl;

    for (vector<int>::iterator ii=v_sol.begin(); ii !=v_sol.end(); ii++) { 
        *ii=*ii+1;  
    }
    Tprint(&v_sol,outputfile);

    cout << " ...and saying goodbye \n";
    Tclean_mat(&AA);  Tclean_mat(&nozeros_AA);
}
