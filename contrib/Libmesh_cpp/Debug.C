#define _IN_DEBUG_DEF_

#include "Debug.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <libgen.h>
#include "libmesh/libmesh.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"

/**
*    Constructor of a Debug object
*    @param fileNamePrefix prefix of the logfile
*    @param comm           MPI Communicator in a multi-process run
*                          (default value MPI_COMM_NULL for a
*                          sequential run)
*/
Debug::Debug() : std::ofstream()
{
}

void Debug::open(const char *out, MPI_Comm comm) 
{
  std::ostringstream s;

  int iproc;
  if (comm == MPI_COMM_NULL)
    s << out << ".log";
  else {
    int iProc;
    MPI_Comm_rank(comm, &iProc);
    s << out << iProc << ".log";
  }
  std::ofstream::open(s.str().c_str());
}

Debug::~Debug()
{
}

Debug & Debug::operator<<(const libMesh::Node &N)
{
  std::ofstream & out = *this;
  out << " id(pid) :" << std::setw(4) << N.id()
	 << "(" << N.processor_id() << ")  xyz = ";
  for (int k=0; k<LIBMESH_DIM; k++)
    out << " " << std::setw(12) << N(k);
  out << '\n';
  out.flush();

  return *this;
}

Debug & Debug::operator<<(const libMesh::Elem &E)
{
  std::ofstream & out = *this;

  out << "id(pid) :" << std::setw(4) << E.id() 
      << "(" << E.processor_id() << ")";
  out << " nodes = ";
  for (unsigned int k=0; k<E.n_nodes(); k++)
    out << " " << std::setw(4) << E.get_node(k)->id()
	<< "(" << E.get_node(k)->processor_id() << ")";
  
  out << " neighbors = ";
  for (unsigned int k=0; k<E.n_neighbors(); k++) {
    libMesh::Elem * EN = E.neighbor(k);
    if (EN)
      if (EN->is_remote())
	out << " remote";
      else
	out << std::setw(4) << EN->id() << "(" << EN->processor_id() << ")";
    else
      out << "   NULL";
  }
  out << " n_systems = " << E.n_systems() << '\n';
  return *this;
}

void Debug::position(const char *file, int line)
{
  std::ofstream & out = *this;
  char * s = strdup(file);
  out << std::right << std::setw(25) 
      << basename(s) << '[' << line << "] ";
  free(s);
}
