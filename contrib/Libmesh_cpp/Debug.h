/*!
*     @file Debug.h
*     @brief definition of the Debug class, used to write text output
*     in a different file for each process in a multi-processes run
*
*     @author Marc Tajchman (marc.tajchman@cea.fr)
*
*     @date 2/12/2011
*/

#ifndef __DEBUG_HXX__
#define __DEBUG_HXX__

#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include "MEDCouplingMemArray.hxx"

namespace libMesh {
  class Node;
  class Elem;
}

class Debug : public std::ofstream {

public:

  Debug();
  void open(const char *fileNamePrefix,  MPI_Comm comm = MPI_COMM_NULL);
  ~Debug();

  void position(const char *, int);
  Debug & operator<< (const libMesh::Node &);
  Debug & operator<< (const libMesh::Elem &);

  typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
  typedef CoutType& (*StandardEndLine)(CoutType&);

  Debug& operator<<(StandardEndLine manip) {
    std::ofstream & f = *this;
    manip(f);
    f.flush();
    return *this;
  }

  template <typename T>
  void array(T *d, long n, long m, const char *s) {
    long i,j;
    std::ofstream & f = *this;
    f << s << '\n';
    for (i=0; i<n; i++) {
      f << " " << std::setw(6) << i << ":";
      for (j=0; j<m; j++)
	f  << " " << std::setw(14) << *d++;
      f << '\n';
    }
  }

  template <typename T>
  void array(std::vector<T *> & d)
  {
    unsigned int j;    
    for (j=0; j<d.size(); j++)
      *this << *(d[j])++ << std::endl;
  }

  template <typename T, typename U>
  void array(const std::set<T *, U> & d)
  {
    typename std::set<T *, U>::const_iterator it = d.begin();
    for (; it != d.end(); it++)
      *this << *(*it);
    *this << std::endl;
  }

  template <typename T>
  void array(std::vector<T> *d, long n, long m, const char *s) {
    long i,j;
    std::ofstream & f = *this;
    f << s << '\n';
    for (i=0; i<n; i++) {
      f << " " << std::setw(6) << i << ":";
      for (j=0; j<m; j++)
	f  << " " << std::setw(14) << *d++;
      f << '\n';
    }
  }

  template <typename T>
  void array(std::vector<std::vector<T> * > & d, 
	     int n, std::vector<int> & m, const char *s) {
    long i,j,p;
    std::ofstream & f = *this;
    f << s << '\n';
    for (i=0; i<n; i++) {
      f << " " << std::setw(6) << i << ":";
      for (p=0; p<d.size(); p++) {
        f << " ";
        for (j=i*m[p]; j<(i+1)*m[p]; j++)
	  f  << " " << std::setw(14) << (*d[p])[j];
      }
      f << '\n';
    }
  }

  void array(const ParaMEDMEM::DataArrayDouble * d, 
	     const char *s) {
    long i,j,p;
    std::ofstream & f = *this;
    f << s << '\n';
    for (i=0; i<d->getNumberOfTuples(); i++) {
      f << " " << std::setw(6) << i << ": ";
      for (j=0; j<d->getNumberOfComponents(); j++)
	f  << " " << std::setw(14) << d->getIJ(i,j);
      f << '\n';
    }
  }

  void array(const ParaMEDMEM::DataArrayInt * d, 
	     const char *s) {
    long i,j,p;
    std::ofstream & f = *this;
    f << s << '\n';
    for (i=0; i<d->getNumberOfTuples(); i++) {
      f << " " << std::setw(6) << i << ": ";
      for (j=0; j<d->getNumberOfComponents(); j++)
	f  << " " << std::setw(14) << d->getIJ(i,j);
      f << '\n';
    }
  }

protected:
  bool _prefix;
};

#ifndef _IN_DEBUG_DEF_
extern
#endif
Debug  fDebug;

#define fDebugPos fDebug.position(__FILE__, __LINE__); fDebug

//template <typename T>
//void arrayDebug(T *d, long n, long m, const char *s)
//{
//  long i,j;
//  fDebug << s << std::endl;
//  for (i=0; i<n; i++) {
//    fDebug << " " << std::setw(6) << i << ":";
//    for (j=0; j<m; j++)
//      fDebug  << " " << std::setw(14) << *d++;
//    fDebug << "\n";
//  }
//}
//
//template <typename T>
//void arrayDebug(std::vector<T *> & d, int n, std::vector<int>& m, const char *s)
//{
//  long i,k;
//  unsigned int j;
//  fDebug << s << std::endl;
//
//  std::vector<T *> e = d;
//  for (i=0; i<n; i++) {
//    fDebug << " " << std::setw(6) << i << ":";
//    for (j=0; j<d.size(); j++)
//      for (k=0; k<m[j]; k++)
//	fDebug << " " << std::setw(14) << *(e[j])++;
//    //    fDebug << std::endl;
//  }
//}
//
#endif

