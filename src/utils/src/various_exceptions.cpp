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


#include <exception>
using namespace std;

#include "various_exceptions.h"


toolsexception::toolsexception() throw() { }

toolsexception::~toolsexception() throw() { }

const char* toolsexception::what() const throw() {
  static const char* what_string = "Tools exception  ";
  return what_string;
  }

const char*  log_negative::what() const throw() {
 static const char* what_string = "Logarithm of a non-positive value  ";
 return what_string;
}


const char* notaprobvector::what() const throw() {
 static const char* what_string = "The vector you are using is not a probability vector";
 return what_string;
}

const char* notamatrix::what() const throw() {
 static const char* what_string = "The argument you are passing is not a matrix";
 return what_string;
}

const char* notajointmatrix::what() const throw() {
 static const char* what_string = "The argument you are passing is not a joint matrix";
 return what_string;
}


const char* dimismatch::what() const throw() {
 static const char* what_string = "Dimension mismatch in the matrices";
 return what_string;
}


const char* emptyrowcolumn::what() const throw() {
 static const char* what_string = "There is an empty row or an empty column";
 return what_string;
}


const char* somethingwronginclustering::what() const throw() {
 static const char* what_string = "Something went wrong in clustering - especially on the distance matrix";
 return what_string;
}

const char* somethingwrongininputmatrix::what() const throw() {
 static const char* what_string = "Something went wrong in input matrix - check the input";
 return what_string;
}

const char* errorduringmerge::what() const throw() {
 static const char* what_string = "an error occured during the merging";
 return what_string;
}
