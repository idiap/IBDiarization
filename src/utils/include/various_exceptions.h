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


#ifndef MATH_EXCEPTION_H 
#define MATH_EXCEPTION_H

#include <exception>

using namespace std;

class toolsexception: public std::exception
{

  public:
      toolsexception() throw();
      ~toolsexception() throw();
      virtual const char* what() const throw();

};



class log_negative: public toolsexception
{
  public:
      const char* what() const throw();
};


class notaprobvector: public toolsexception
{
  public:
      const char* what() const throw();
};

class notamatrix: public toolsexception
{
  public:
      const char* what() const throw();
};

class dimismatch: public toolsexception
{
  public:
      const char* what() const throw();
};

class emptyrowcolumn: public toolsexception
{
  public:
      const char* what() const throw();
};

class notajointmatrix: public toolsexception
{
  public:
      const char* what() const throw();
};

class somethingwronginclustering: public toolsexception
{
  public:
      const char* what() const throw();
};


class somethingwrongininputmatrix: public toolsexception
{
  public:
      const char* what() const throw();
};

class errorduringmerge: public toolsexception
{
  public:
      const char* what() const throw();
};


#endif

