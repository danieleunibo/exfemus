#define FDEBUG 0

#include <iostream>
#include <iomanip>
#include <cstring>

// Libmesh includes
#include "libmesh/point.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_data.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_quad4.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex20.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_quad9.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/remote_elem.h"
#include "libmesh/boundary_info.h"
 
#include "Debug.hxx"
#include "ParallelMeshExtended.hxx"
#include "Timer.hxx"
 
using namespace ParaMEDMEM;
using namespace libMesh;

class TestEqual {
public:
  TestEqual(double eps, int dim) { _eps = eps; _dim = dim; }
  bool operator() (double *a, double *b) {
    bool test = true;
    for (int i=0; i<_dim; i++)
      if (std::fabs(*a++ - *b++) > _eps) { test = false; break; }
    return test;
  }
protected:
  double _eps;
  int _dim;
};

ParallelMeshExtended::ParallelMeshExtended
(const ParaMEDMEM::MEDCouplingUMesh & mesh, 
 MPI_Comm comm)
  : libMesh::ParallelMesh(mesh.getSpaceDimension()),
    _mesh(mesh),
    _comm(comm)
{

  Timer T;

  int iProc, jProc, kProc, nProcs;
  if (_comm != MPI_COMM_NULL) {
    MPI_Comm_size(_comm, &nProcs);
    MPI_Comm_rank(_comm, &iProc);
  }
  else {
    nProcs = 1;
    iProc = 0;
  }

  T.start();

  int i, j, dim = mesh.getSpaceDimension();
  unsigned int ui;
  _meshDim = mesh.getMeshDimension();

  const ParaMEDMEM::DataArrayDouble * LocalNodes = mesh.getCoords();
  int nLocalNodes = LocalNodes->getNumberOfTuples();
  int nLocalElems = mesh.getNumberOfCells();
  int nOwnedNodes = nLocalNodes;

  if (FDEBUG) {
    char t[200];
    sprintf(t, "nodes on %d", iProc);
    fDebug.array(LocalNodes->getConstPointer(), 
		 LocalNodes->getNumberOfTuples(), dim, t);
  }

  _proc_id.resize(nLocalNodes);
  _node_id.resize(nLocalNodes);
  _elem_id.resize(nLocalElems);

  for (i=0; i<nLocalNodes; i++) {
    _proc_id[i] = iProc;
    _node_id[i] = i;
  }
  
  for (i=0; i<nLocalElems; i++) {
    _elem_id[i] = i;
  }

  std::vector< std::vector<int> > index(nProcs);
  std::vector< std::vector<int> > index_id(nProcs);

  std::vector< std::vector<double> > coord(nProcs);

  ParaMEDMEM::DataArrayInt * d = mesh.findBoundaryNodes();
  index[iProc].resize(d->getNbOfElems());
  std::memcpy(&(index[iProc][0]), d->getConstPointer(), d->getNbOfElems());
  d->decrRef();

  for (unsigned in=0; in<index[iProc].size();in++)
    mesh.getCoordinatesOfNode(index[iProc][in], coord[iProc]);
  index_id[iProc].resize(index[iProc].size());

  int *sOut = new int[nProcs];
  int sIn = index[iProc].size();

  if (_comm != MPI_COMM_NULL) {
    MPI_Allgather(&sIn,1,MPI_INT,sOut,1,MPI_INT,_comm);
  }
  for (i=0; i<nProcs; i++) 
    if (i != iProc) {
      index[i].resize(sOut[i]);
      index_id[i].resize(sOut[i]);
      coord[i].resize(sOut[i] * dim);
    }

  if (_comm != MPI_COMM_NULL) {
    for (i=0; i<nProcs; i++) {
      MPI_Bcast(&(index[i][0]), index[i].size(), MPI_INT, i, _comm);
      MPI_Bcast(&(coord[i][0]), coord[i].size(), MPI_DOUBLE, i, _comm);
    }

    if (FDEBUG) {
      fDebug << "index[" << i << "] =";
      for (unsigned int j=0; j<index[i].size(); j++)
	fDebug << " " << index[i][j];
      fDebug << std::endl;
    }
  }

  fDebug << "Elapsed     (exchange boundaries) " << T.elapsed() << std::endl;

  T.start();
  if (iProc > 0) {
    
    int i,j;
    TestEqual T(1e-12, dim);
    double * c_iProc = &(coord[iProc][0]);
    int n_iProc = sOut[iProc];

    for (jProc = iProc-1; jProc >=0; jProc--) {
      double * c_jProc = &(coord[jProc][0]);
      int n_jProc = sOut[jProc];

      for (i=0; i<n_iProc; i++) {

	bool Found = false;
	double * a_i = c_iProc + i*dim;
	for (j=0; j<n_jProc; j++) {
	  double * b_j = c_jProc + j*dim;
	  if (T(a_i, b_j)) {

	    if (FDEBUG) {
	      int k;
	      fDebug << "a_" << i;
	      for (k=0; k<dim; k++)
		fDebug << " " << a_i[k];
	      fDebug << std::endl;
	      fDebug << "b_" << j;
	      for (k=0; k<dim; k++)
		fDebug << " " << b_j[k];
	      fDebug << std::endl;
	    }

	    Found = true; break;
	  }
	}
	if (Found) {
	  int id = index[iProc][i];
	  _node_id[id] = index[jProc][j];
	  _proc_id[id] = jProc;
	}
      }
      
    }
  }

  nOwnedNodes = 0;
  for (i=0; i<nLocalNodes; i++)
    if (_proc_id[i] == iProc) nOwnedNodes++;

  if (FDEBUG) {
    std::vector<std::vector<int> *> v(2); v[0] = &_node_id ; v[1] = &_proc_id;
    std::vector<int>   m(2); m[0] = 1; m[1] = 1;
    fDebug.array(v, nLocalNodes, m, "Node id - proc :");
  }

  int * nodeIdOffset = new int[nProcs+1];
  int * elemIdOffset = new int[nProcs+1];
  for (i=0; i<nProcs+1; i++)
    nodeIdOffset[i] = elemIdOffset[i] = 0;

  if (_comm != MPI_COMM_NULL) {
    MPI_Allgather(&nLocalElems,      1, MPI_INT,elemIdOffset+1, 1, 
		  MPI_INT, _comm);
    MPI_Allgather(&nOwnedNodes, 1, MPI_INT,nodeIdOffset+1, 1, 
		  MPI_INT, _comm);
  }
  if (FDEBUG) {
    fDebug.array(nodeIdOffset, 1, nProcs+1, "Node offsets :");
    fDebug.array(elemIdOffset, 1, nProcs+1, "Element offsets :");
    fDebug << "Elapsed     (node offsets)        " << T.elapsed() << std::endl;
  }

  T.start();

  for (i=1; i<=nProcs; i++) {
    nodeIdOffset[i] += nodeIdOffset[i-1];
    elemIdOffset[i] += elemIdOffset[i-1];
  }
  
  if (FDEBUG) {
    fDebug.array(nodeIdOffset, 1, nProcs+1, "Node offsets :");
    fDebug.array(elemIdOffset, 1, nProcs+1, "Element offsets :");
  }
  
  int iNode = 0;
  for (int i=0; i<nLocalNodes; i++) {
    if (_proc_id[i] == iProc)
      _node_id[i] = nodeIdOffset[_proc_id[i]] + iNode++;
  }

  for (int i=0; i<nLocalElems; i++) {
    _elem_id[i] += elemIdOffset[iProc];
  }

  for (unsigned int k=0; k<index[iProc].size(); k++)
    index_id[iProc][k] = _proc_id[index[iProc][k]] == iProc 
      ? _node_id[index[iProc][k]] 
      : -1;
 
  if (_comm != MPI_COMM_NULL) {
    MPI_Barrier(_comm);
    for (jProc=0; jProc<nProcs; jProc++) {
      MPI_Bcast(&(index_id[jProc][0]), index_id[jProc].size(), 
		MPI_INT, jProc, _comm);
    }
    MPI_Barrier(_comm);
  }
  fDebug << "Elapsed     (elements offsets)    " << T.elapsed() << std::endl;
 
  if (FDEBUG) 
    {
      fDebug << "ok ici " << std::endl;
      for (jProc=0; jProc<nProcs; jProc++) {
	std::vector<std::vector<int> *> v(2); 
	v[0] = &(index[jProc]); 
	v[1] = &(index_id[jProc]);
	std::vector<int>   m(2); m[0] = 1; m[1] = 1;
	char s[200];
	sprintf(s, "Index - Index_id %d", jProc);
	fDebug.array(v, index[jProc].size(), m, s);
      }
    }

  for (unsigned int i=0; i<index[iProc].size(); i++) {
    j    = index[iProc][i];
    jProc = _proc_id[j];
    if (FDEBUG) 
      fDebug << "j = " << j
	     << " jProc = " << jProc 
	     << std::endl;

    if (jProc != iProc) {
      int j_id = _node_id[j];
      if (FDEBUG) 
	fDebug << std::endl 
	       << " i = " << i 
	       << " j_id = " << j_id 
	       << " jProc = " << jProc 
	       << std::endl;
      for (unsigned int k=0; k<index[jProc].size(); k++)
	if (index[jProc][k] == j_id) 
	  _node_id[j] = index_id[jProc][k];
      if (FDEBUG)
	fDebug << "_node_id[j] = " << _node_id[j] << std::endl;
    }
  }

  if (FDEBUG) {
    fDebug << std::endl;
    std::vector<std::vector<int> *> v(2); 
    v[0] = &_node_id; 
    v[1] = &_proc_id;
    std::vector<int>   m(2); m[0] = 1; m[1] = 1;
    fDebug.array(v, nLocalNodes, m, "Node id - proc :");
  }

  Point xyz;
  int id;
  const double * c = LocalNodes->getConstPointer();
  for (int i=0; i<nLocalNodes; i++) {
    xyz(0) = *c++;
    xyz(1) = (dim > 1) ? *c++ : 0.0;
    xyz(2) = (dim > 2) ? *c++ : 0.0;
    id = _node_id[i];
    this->add_point (xyz, id, _proc_id[i]);

    if (FDEBUG)
      fDebug << "point (" 
	     << std::setw(9) << xyz(0) << ", " 
	     << std::setw(9) << xyz(1) << ", " 
	     << std::setw(9) << xyz(2) << ")"
	     << ", id = " << std::setw(6) << id 
	     << " _proc_id " << std::setw(6) << _proc_id[i] << std::endl;
  }  

  const DataArrayInt * meshConnect = mesh.getNodalConnectivity();
  const DataArrayInt * meshConnectIdx = mesh.getNodalConnectivityIndex();
  
  for (int iElem=0; iElem < nLocalElems; iElem++) {

    INTERP_KERNEL::NormalizedCellType t 
      = mesh.getTypeOfCell (iElem);
    Elem * elem;
    Elem * elem1;
     
  
    int n_nodes_elem;
    int n_nodes_elem1;
    

    switch (t) {
    case INTERP_KERNEL::NORM_SEG2:
      elem = new Edge2;
      n_nodes_elem = 2;
      break;
    case INTERP_KERNEL::NORM_TRI3:
      elem = new Tri3;
      n_nodes_elem = 3;
      break;
    case INTERP_KERNEL::NORM_TETRA4:
      elem = new Tet4;
      n_nodes_elem = 4;
      break;
    case INTERP_KERNEL::NORM_TETRA10:
      elem = new Tet10;
      n_nodes_elem = 10;
      break; 
    case INTERP_KERNEL::NORM_HEXA8:
       elem = new Hex8;
       n_nodes_elem = 8;
        break;
    case INTERP_KERNEL::NORM_HEXA20:
      elem = new Hex20;
      n_nodes_elem = 20;
      break;   
    case INTERP_KERNEL::NORM_HEXA27:
      elem = new Hex27;
      n_nodes_elem = 27;
      break; 
    case INTERP_KERNEL::NORM_QUAD4:
      elem = new Quad4;
      n_nodes_elem = 4;
      break;
    case INTERP_KERNEL::NORM_QUAD8:
      elem = new Quad8;
      n_nodes_elem = 8;
      break;
       case INTERP_KERNEL::NORM_QUAD9:
      elem = new Quad9;
      n_nodes_elem = 9;
      break;
    default:
      elem = NULL;
      n_nodes_elem = 0;
    }
    
    
    int cstart = meshConnectIdx->getIJ(iElem, 0);
    //   int cend = meshConnectIdx->getIJ(iElem+1, 0);

    int transl[20]={0,3,2,1,4,7,6,5,11,10,9,8,19,18,17,16,12,15,14,13};
    int transl27[27]={0,3,2,1,4,7,6,5,11,10,9,8,19,18,17,16,12,15,14,13,20,24,23,22,21,25,26};
    int transl9[9]={1,0,3,2,4,7,6,5,8};
    
    for (int i=0; i<n_nodes_elem; i++)
      elem->set_node(i) 
	= node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
 
        
        
        
      
    switch (t) {
    case INTERP_KERNEL::NORM_SEG2:
      if (elem->volume() < 0.0) {
        Node * n = elem->get_node(0);
        elem->set_node(0) = elem->get_node(1);
        elem->set_node(1) = n;
      }
      break;
    case INTERP_KERNEL::NORM_TRI3:
    case INTERP_KERNEL::NORM_QUAD4:
      if (elem->volume() < 0.0) {
        Node * n;
        for (int i=0; i<n_nodes_elem/2; i++) {
          n = elem->get_node(i);
          elem->set_node(i) = elem->get_node(n_nodes_elem-1-i);
          elem->set_node(n_nodes_elem-1-i) = n;
        }
      }
      break;
      case INTERP_KERNEL::NORM_QUAD8:
      if (elem->volume() < 0.0) {
        Node * n;
        for (int i=0; i<2; i++) {
          n = elem->get_node(i);
          elem->set_node(i) = elem->get_node(4-1-i);
          elem->set_node(4-1-i) = n;
        }
         n = elem->get_node(5);
         elem->set_node(5) = elem->get_node(7);
          elem->set_node(7) = n;
         
      } 
      break;
      case INTERP_KERNEL::NORM_QUAD9:
        
    for (int i=0; i<n_nodes_elem; i++)
      elem->set_node(transl9[i]) 
        = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
        
      if (elem->volume() < 0.0) {
        Node * n;
        for (int i=0; i<2; i++) {
          n = elem->get_node(i);
          elem->set_node(i) = elem->get_node(4-1-i);
          elem->set_node(4-1-i) = n;
        }
         n = elem->get_node(5);
         elem->set_node(5) = elem->get_node(7);
          elem->set_node(7) = n;
         
      } 
      
      
      break; 
    case INTERP_KERNEL::NORM_TETRA4:
      if (elem->volume() < 0.0) {
        Node * n = elem->get_node(1);
        elem->set_node(1) = elem->get_node(2);
        elem->set_node(2) = n;
      }
     break; 
     case INTERP_KERNEL::NORM_TETRA10:
      
      elem1 = new Tet4;
      n_nodes_elem1 = 4;
        for (int i=0; i<n_nodes_elem1; i++)
      elem1->set_node(i) = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]); 
       
      if (elem1->volume() < 0.0) {
        Node * n = elem->get_node(1);
        elem->set_node(1) = elem->get_node(2);
        elem->set_node(2) = n;
        n = elem->get_node(4);
        elem->set_node(4) = elem->get_node(6);
        elem->set_node(6) = n;
        n = elem->get_node(9);
        elem->set_node(9) = elem->get_node(8);
        elem->set_node(8) = n;
       } 
       break;
     case INTERP_KERNEL::NORM_HEXA8:
       
       if (elem->volume() < 0.0) {
        Node * n = elem->get_node(1);
        elem->set_node(1) = elem->get_node(3);
        elem->set_node(3) = n;
        n = elem->get_node(5);
        elem->set_node(5) = elem->get_node(7);
        elem->set_node(7) = n;

      }
      break;
       
     case INTERP_KERNEL::NORM_HEXA20:
       
       
       for (int i=0; i<n_nodes_elem; i++)
        elem->set_node(transl[i]) 
        = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
       
       
       elem1 = new Hex8;
       n_nodes_elem1 = 8;
       for (int i=0; i<n_nodes_elem1; i++)
      elem1->set_node(transl[i]) = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]); 
      
      
      
      
       if (elem1->volume() < 0.0) {
         
         
        Node * n = elem->get_node(1);
        elem->set_node(1) = elem->get_node(3);
        elem->set_node(3) = n;
        n = elem->get_node(8);
        elem->set_node(8) = elem->get_node(11);
        elem->set_node(11) = n;
        n = elem->get_node(9);
        elem->set_node(9) = elem->get_node(10);
        elem->set_node(10) = n;
        n = elem->get_node(5);
        elem->set_node(5) = elem->get_node(7);
        elem->set_node(7) = n;
        n = elem->get_node(12);
        elem->set_node(12) = elem->get_node(19);
        elem->set_node(19) = n;
        n = elem->get_node(13);
        elem->set_node(13) = elem->get_node(18);
        elem->set_node(18) = n;
        n = elem->get_node(14);
        elem->set_node(14) = elem->get_node(17);
        elem->set_node(17) = n;
        
      }
      
      break;
      
       case INTERP_KERNEL::NORM_HEXA27:
       
       
       for (int i=0; i<n_nodes_elem; i++)
        elem->set_node(transl27[i]) 
        = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
       
       
       
       elem1= new Hex8;
       n_nodes_elem1 = 8;
       for (int i=0; i<n_nodes_elem1; i++)
      elem1->set_node(transl27[i]) = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);       

       if (elem1->volume() < 0.0) {
        Node * n = elem->get_node(1);
        elem->set_node(1) = elem->get_node(3);
        elem->set_node(3) = n;
        n = elem->get_node(8);
        elem->set_node(8) = elem->get_node(11);
        elem->set_node(11) = n;
        n = elem->get_node(9);
        elem->set_node(9) = elem->get_node(10);
        elem->set_node(10) = n;
        n = elem->get_node(5);
        elem->set_node(5) = elem->get_node(7);
        elem->set_node(7) = n;
        n = elem->get_node(12);
        elem->set_node(12) = elem->get_node(19);
        elem->set_node(19) = n;
        n = elem->get_node(13);
        elem->set_node(13) = elem->get_node(18);
        elem->set_node(18) = n;
        n = elem->get_node(14);
        elem->set_node(14) = elem->get_node(17);
        elem->set_node(17) = n;
        n = elem->get_node(22);
        elem->set_node(22) = elem->get_node(23);
        elem->set_node(23) = n;
        n = elem->get_node(24);
        elem->set_node(24) = elem->get_node(21);
        elem->set_node(21) = n;
        
      }
     
      break;
    default:
      break;
    }  

   
    elem->set_id(_elem_id[iElem]);
    
    elem->processor_id() = iProc;

    this->add_elem (elem);

    if (FDEBUG) {
      fDebug << "element (";
      for (int i=0; i<n_nodes_elem; i++)
	fDebug << std::setw(6) << elem->get_node(i)->id() << " ";
      fDebug     << "), id = " << std::setw(6) << _elem_id[iElem]
		 << " _proc_id " << std::setw(6) << iProc << std::endl;
    }

  }

  delete [] sOut;
  fDebug << "Elapsed     (insert nodes-elems)    " << T.elapsed() << std::endl;

  _n_nodes = nodeIdOffset[nProcs];
  _n_elem = elemIdOffset[nProcs];
  if (FDEBUG) {
    fDebug << "n_nodes : " << _n_nodes << " " << this->n_local_nodes() << std::endl;
    fDebug << "n_nodes : " << this->n_nodes() << " " << this->n_local_nodes() << std::endl;
  }

  _is_serial = false;
  _partitioner = AutoPtr<libMesh::Partitioner>(NULL);
  T.start();
  //  MeshCommunication().gather_neighboring_elements(*this);
  prepare_for_use(false);
  fDebug << "Elapsed     (prepare)    " << T.elapsed() << std::endl;

  //  fDebug << boundary_info->n_boundary_ids() << std::endl;

  int nBdyFaces = 0;
  MeshBase::const_element_iterator       el     
    = active_local_elements_begin();
  const MeshBase::const_element_iterator end_el 
    = active_local_elements_end();
  for ( ; el != end_el; ++el) {
    const Elem* elem = *el;
    for (unsigned int side=0; side<elem->n_sides(); side++)
      if (elem->neighbor(side) == NULL)
	{
	  nBdyFaces++;
	  if (FDEBUG) fDebug << " elem " << elem->id() 
		 << " side : " << side << std::endl;
	}
  }

  delete [] nodeIdOffset;
  delete [] elemIdOffset;
}

ParallelMeshExtended::~ParallelMeshExtended()
{
}
