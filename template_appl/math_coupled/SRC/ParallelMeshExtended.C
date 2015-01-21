#define FDEBUG 0

#include <iostream>
#include <iomanip>
#include <cstring>

// Libmesh includes
// #include "libmesh/point.h"
// #include "libmesh/mesh.h"
// #include "libmesh/mesh_base.h"
// #include "libmesh/mesh_data.h"
// #include "libmesh/parallel.h"
// #include "libmesh/parallel_mesh.h"
// #include "libmesh/cell_tet4.h"
// #include "libmesh/cell_tet10.h"
// #include "libmesh/face_tri3.h"
// #include "libmesh/face_quad4.h"
// #include "libmesh/cell_hex8.h"
// #include "libmesh/cell_hex20.h"
// #include "libmesh/cell_hex27.h"
// #include "libmesh/face_quad8.h"
// #include "libmesh/face_quad9.h"
// #include "libmesh/edge_edge2.h"
// #include "libmesh/mesh_communication.h"
// #include "libmesh/remote_elem.h"
// #include "libmesh/boundary_info.h"

// #include "Debug.hxx"
#include "ParallelMeshExtended.h"
#include "MGGeomEl.h"
// #include "Timer.hxx"
// #ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDFileMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
// #endif
#include "MGUtils.h"
// using namespace ParaMEDMEM;
// using namespace libMesh;

class TestEqual {
public:
  TestEqual(double eps, int dim) { _eps = eps; _dim = dim; }
  bool operator()(double *a, double *b) {
    bool test = true;
    for(int i=0; i<_dim; i++)
      if(std::fabs(*a++ - *b++) > _eps) { test = false; break; }
    return test;
  }
protected:
  double _eps;
  int _dim;
};

// ===============================================================
// This is the constructor
ParallelMeshExtended::ParallelMeshExtended(
  const ParaMEDMEM::MEDCouplingUMesh & mesh,
  const ParallelObjectM              & comm,
  MGUtils                            & mgutils,
  MGGeomEl                           & geomel,
  MPI_Comm                           comm1
) : MGMesh(
    comm,
    mgutils,
    geomel
  ),
  _num_meshes(1),
  _mesh(mesh),
  _comm(comm1) {

//   Timer T;
//
  int iProc, jProc, kProc, nProcs;
  if(_comm != MPI_COMM_NULL) {
    MPI_Comm_size(_comm, &nProcs);
    MPI_Comm_rank(_comm, &iProc);
  } else {
    nProcs = 1;
    iProc = 0;
  }
//   _meshes.resize(_num_meshes);
 read_bc_id(_NoLevels-1);
  read_mat_id(_NoLevels-1);
//    _meshes[0]=&mesh;
//
//   T.start();
//
  int i, j, dim = mesh.getSpaceDimension();
  unsigned int ui;
  _meshDim = mesh.getMeshDimension();
//
//   const ParaMEDMEM::DataArrayDouble * LocalNodes = mesh.getCoords();
//   int nLocalNodes = LocalNodes->getNumberOfTuples();
//   int nLocalElems = mesh.getNumberOfCells();
//   int nOwnedNodes = nLocalNodes;
//
//   if (FDEBUG) {
//     char t[200];
//     sprintf(t, "nodes on %d", iProc);
//     fDebug.array(LocalNodes->getConstPointer(),
// 		 LocalNodes->getNumberOfTuples(), dim, t);
//   }
//
//   _proc_id.resize(nLocalNodes);
//   _node_id.resize(nLocalNodes);
//   _elem_id.resize(nLocalElems);
//
//   for (i=0; i<nLocalNodes; i++) {
//     _proc_id[i] = iProc;
//     _node_id[i] = i;
//   }
//
//   for (i=0; i<nLocalElems; i++) {
//     _elem_id[i] = i;
//   }
//
//   std::vector< std::vector<int> > index(nProcs);
//   std::vector< std::vector<int> > index_id(nProcs);
//
//   std::vector< std::vector<double> > coord(nProcs);
//
//   ParaMEDMEM::DataArrayInt * d = mesh.findBoundaryNodes();
//   index[iProc].resize(d->getNbOfElems());
//   std::memcpy(&(index[iProc][0]), d->getConstPointer(), d->getNbOfElems());
//   d->decrRef();
//
//   for (unsigned in=0; in<index[iProc].size();in++)
//     mesh.getCoordinatesOfNode(index[iProc][in], coord[iProc]);
//   index_id[iProc].resize(index[iProc].size());
//
//   int *sOut = new int[nProcs];
//   int sIn = index[iProc].size();
//
//   if (_comm != MPI_COMM_NULL) {
//     MPI_Allgather(&sIn,1,MPI_INT,sOut,1,MPI_INT,_comm);
//   }
//   for (i=0; i<nProcs; i++)
//     if (i != iProc) {
//       index[i].resize(sOut[i]);
//       index_id[i].resize(sOut[i]);
//       coord[i].resize(sOut[i] * dim);
//     }
//
//   if (_comm != MPI_COMM_NULL) {
//     for (i=0; i<nProcs; i++) {
//       MPI_Bcast(&(index[i][0]), index[i].size(), MPI_INT, i, _comm);
//       MPI_Bcast(&(coord[i][0]), coord[i].size(), MPI_DOUBLE, i, _comm);
//     }
//
//     if (FDEBUG) {
//       fDebug << "index[" << i << "] =";
//       for (unsigned int j=0; j<index[i].size(); j++)
// 	fDebug << " " << index[i][j];
//       fDebug << std::endl;
//     }
//   }
//
//   fDebug << "Elapsed     (exchange boundaries) " << T.elapsed() << std::endl;
//
//   T.start();
//   if (iProc > 0) {
//
//     int i,j;
//     TestEqual T(1e-12, dim);
//     double * c_iProc = &(coord[iProc][0]);
//     int n_iProc = sOut[iProc];
//
//     for (jProc = iProc-1; jProc >=0; jProc--) {
//       double * c_jProc = &(coord[jProc][0]);
//       int n_jProc = sOut[jProc];
//
//       for (i=0; i<n_iProc; i++) {
//
// 	bool Found = false;
// 	double * a_i = c_iProc + i*dim;
// 	for (j=0; j<n_jProc; j++) {
// 	  double * b_j = c_jProc + j*dim;
// 	  if (T(a_i, b_j)) {
//
// 	    if (FDEBUG) {
// 	      int k;
// 	      fDebug << "a_" << i;
// 	      for (k=0; k<dim; k++)
// 		fDebug << " " << a_i[k];
// 	      fDebug << std::endl;
// 	      fDebug << "b_" << j;
// 	      for (k=0; k<dim; k++)
// 		fDebug << " " << b_j[k];
// 	      fDebug << std::endl;
// 	    }
//
// 	    Found = true; break;
// 	  }
// 	}
// 	if (Found) {
// 	  int id = index[iProc][i];
// 	  _node_id[id] = index[jProc][j];
// 	  _proc_id[id] = jProc;
// 	}
//       }
//
//     }
//   }
//
//   nOwnedNodes = 0;
//   for (i=0; i<nLocalNodes; i++)
//     if (_proc_id[i] == iProc) nOwnedNodes++;
//
//   if (FDEBUG) {
//     std::vector<std::vector<int> *> v(2); v[0] = &_node_id ; v[1] = &_proc_id;
//     std::vector<int>   m(2); m[0] = 1; m[1] = 1;
//     fDebug.array(v, nLocalNodes, m, "Node id - proc :");
//   }
//
//   int * nodeIdOffset = new int[nProcs+1];
//   int * elemIdOffset = new int[nProcs+1];
//   for (i=0; i<nProcs+1; i++)
//     nodeIdOffset[i] = elemIdOffset[i] = 0;
//
//   if (_comm != MPI_COMM_NULL) {
//     MPI_Allgather(&nLocalElems,      1, MPI_INT,elemIdOffset+1, 1,
// 		  MPI_INT, _comm);
//     MPI_Allgather(&nOwnedNodes, 1, MPI_INT,nodeIdOffset+1, 1,
// 		  MPI_INT, _comm);
//   }
//   if (FDEBUG) {
//     fDebug.array(nodeIdOffset, 1, nProcs+1, "Node offsets :");
//     fDebug.array(elemIdOffset, 1, nProcs+1, "Element offsets :");
//     fDebug << "Elapsed     (node offsets)        " << T.elapsed() << std::endl;
//   }
//
//   T.start();
//
//   for (i=1; i<=nProcs; i++) {
//     nodeIdOffset[i] += nodeIdOffset[i-1];
//     elemIdOffset[i] += elemIdOffset[i-1];
//   }
//
//   if (FDEBUG) {
//     fDebug.array(nodeIdOffset, 1, nProcs+1, "Node offsets :");
//     fDebug.array(elemIdOffset, 1, nProcs+1, "Element offsets :");
//   }
//
//   int iNode = 0;
//   for (int i=0; i<nLocalNodes; i++) {
//     if (_proc_id[i] == iProc)
//       _node_id[i] = nodeIdOffset[_proc_id[i]] + iNode++;
//   }
//
//   for (int i=0; i<nLocalElems; i++) {
//     _elem_id[i] += elemIdOffset[iProc];
//   }
//
//   for (unsigned int k=0; k<index[iProc].size(); k++)
//     index_id[iProc][k] = _proc_id[index[iProc][k]] == iProc
//       ? _node_id[index[iProc][k]]
//       : -1;
//
//   if (_comm != MPI_COMM_NULL) {
//     MPI_Barrier(_comm);
//     for (jProc=0; jProc<nProcs; jProc++) {
//       MPI_Bcast(&(index_id[jProc][0]), index_id[jProc].size(),
// 		MPI_INT, jProc, _comm);
//     }
//     MPI_Barrier(_comm);
//   }
//   fDebug << "Elapsed     (elements offsets)    " << T.elapsed() << std::endl;
//
//   if (FDEBUG)
//     {
//       fDebug << "ok ici " << std::endl;
//       for (jProc=0; jProc<nProcs; jProc++) {
// 	std::vector<std::vector<int> *> v(2);
// 	v[0] = &(index[jProc]);
// 	v[1] = &(index_id[jProc]);
// 	std::vector<int>   m(2); m[0] = 1; m[1] = 1;
// 	char s[200];
// 	sprintf(s, "Index - Index_id %d", jProc);
// 	fDebug.array(v, index[jProc].size(), m, s);
//       }
//     }
//
//   for (unsigned int i=0; i<index[iProc].size(); i++) {
//     j    = index[iProc][i];
//     jProc = _proc_id[j];
//     if (FDEBUG)
//       fDebug << "j = " << j
// 	     << " jProc = " << jProc
// 	     << std::endl;
//
//     if (jProc != iProc) {
//       int j_id = _node_id[j];
//       if (FDEBUG)
// 	fDebug << std::endl
// 	       << " i = " << i
// 	       << " j_id = " << j_id
// 	       << " jProc = " << jProc
// 	       << std::endl;
//       for (unsigned int k=0; k<index[jProc].size(); k++)
// 	if (index[jProc][k] == j_id)
// 	  _node_id[j] = index_id[jProc][k];
//       if (FDEBUG)
// 	fDebug << "_node_id[j] = " << _node_id[j] << std::endl;
//     }
//   }
//
//   if (FDEBUG) {
//     fDebug << std::endl;
//     std::vector<std::vector<int> *> v(2);
//     v[0] = &_node_id;
//     v[1] = &_proc_id;
//     std::vector<int>   m(2); m[0] = 1; m[1] = 1;
//     fDebug.array(v, nLocalNodes, m, "Node id - proc :");
//   }
//
//   Point xyz;
//   int id;
//   const double * c = LocalNodes->getConstPointer();
//   for (int i=0; i<nLocalNodes; i++) {
//     xyz(0) = *c++;
//     xyz(1) = (dim > 1) ? *c++ : 0.0;
//     xyz(2) = (dim > 2) ? *c++ : 0.0;
//     id = _node_id[i];
//     this->add_point (xyz, id, _proc_id[i]);
//
//     if (FDEBUG)
//       fDebug << "point ("
// 	     << std::setw(9) << xyz(0) << ", "
// 	     << std::setw(9) << xyz(1) << ", "
// 	     << std::setw(9) << xyz(2) << ")"
// 	     << ", id = " << std::setw(6) << id
// 	     << " _proc_id " << std::setw(6) << _proc_id[i] << std::endl;
//   }
//
//   const DataArrayInt * meshConnect = mesh.getNodalConnectivity();
//   const DataArrayInt * meshConnectIdx = mesh.getNodalConnectivityIndex();
//
//   for (int iElem=0; iElem < nLocalElems; iElem++) {
//
//     INTERP_KERNEL::NormalizedCellType t
//       = mesh.getTypeOfCell (iElem);
//     Elem * elem;
//     Elem * elem1;
//
//
//     int n_nodes_elem;
//     int n_nodes_elem1;
//
//
//     switch (t) {
//     case INTERP_KERNEL::NORM_SEG2:
//       elem = new Edge2;
//       n_nodes_elem = 2;
//       break;
//     case INTERP_KERNEL::NORM_TRI3:
//       elem = new Tri3;
//       n_nodes_elem = 3;
//       break;
//     case INTERP_KERNEL::NORM_TETRA4:
//       elem = new Tet4;
//       n_nodes_elem = 4;
//       break;
//     case INTERP_KERNEL::NORM_TETRA10:
//       elem = new Tet10;
//       n_nodes_elem = 10;
//       break;
//     case INTERP_KERNEL::NORM_HEXA8:
//        elem = new Hex8;
//        n_nodes_elem = 8;
//         break;
//     case INTERP_KERNEL::NORM_HEXA20:
//       elem = new Hex20;
//       n_nodes_elem = 20;
//       break;
//     case INTERP_KERNEL::NORM_HEXA27:
//       elem = new Hex27;
//       n_nodes_elem = 27;
//       break;
//     case INTERP_KERNEL::NORM_QUAD4:
//       elem = new Quad4;
//       n_nodes_elem = 4;
//       break;
//     case INTERP_KERNEL::NORM_QUAD8:
//       elem = new Quad8;
//       n_nodes_elem = 8;
//       break;
//        case INTERP_KERNEL::NORM_QUAD9:
//       elem = new Quad9;
//       n_nodes_elem = 9;
//       break;
//     default:
//       elem = NULL;
//       n_nodes_elem = 0;
//     }
//
//
//     int cstart = meshConnectIdx->getIJ(iElem, 0);
//     //   int cend = meshConnectIdx->getIJ(iElem+1, 0);
//
//     int transl[20]={0,3,2,1,4,7,6,5,11,10,9,8,19,18,17,16,12,15,14,13};
//     int transl27[27]={0,3,2,1,4,7,6,5,11,10,9,8,19,18,17,16,12,15,14,13,20,24,23,22,21,25,26};
//     int transl9[9]={1,0,3,2,4,7,6,5,8};
//
//     for (int i=0; i<n_nodes_elem; i++)
//       elem->set_node(i)
// 	= node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
//
//
//
//
//
//     switch (t) {
//     case INTERP_KERNEL::NORM_SEG2:
//       if (elem->volume() < 0.0) {
//         Node * n = elem->get_node(0);
//         elem->set_node(0) = elem->get_node(1);
//         elem->set_node(1) = n;
//       }
//       break;
//     case INTERP_KERNEL::NORM_TRI3:
//     case INTERP_KERNEL::NORM_QUAD4:
//       if (elem->volume() < 0.0) {
//         Node * n;
//         for (int i=0; i<n_nodes_elem/2; i++) {
//           n = elem->get_node(i);
//           elem->set_node(i) = elem->get_node(n_nodes_elem-1-i);
//           elem->set_node(n_nodes_elem-1-i) = n;
//         }
//       }
//       break;
//       case INTERP_KERNEL::NORM_QUAD8:
//       if (elem->volume() < 0.0) {
//         Node * n;
//         for (int i=0; i<2; i++) {
//           n = elem->get_node(i);
//           elem->set_node(i) = elem->get_node(4-1-i);
//           elem->set_node(4-1-i) = n;
//         }
//          n = elem->get_node(5);
//          elem->set_node(5) = elem->get_node(7);
//           elem->set_node(7) = n;
//
//       }
//       break;
//       case INTERP_KERNEL::NORM_QUAD9:
//
//     for (int i=0; i<n_nodes_elem; i++)
//       elem->set_node(transl9[i])
//         = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
//
//       if (elem->volume() < 0.0) {
//         Node * n;
//         for (int i=0; i<2; i++) {
//           n = elem->get_node(i);
//           elem->set_node(i) = elem->get_node(4-1-i);
//           elem->set_node(4-1-i) = n;
//         }
//          n = elem->get_node(5);
//          elem->set_node(5) = elem->get_node(7);
//           elem->set_node(7) = n;
//
//       }
//
//
//       break;
//     case INTERP_KERNEL::NORM_TETRA4:
//       if (elem->volume() < 0.0) {
//         Node * n = elem->get_node(1);
//         elem->set_node(1) = elem->get_node(2);
//         elem->set_node(2) = n;
//       }
//      break;
//      case INTERP_KERNEL::NORM_TETRA10:
//
//       elem1 = new Tet4;
//       n_nodes_elem1 = 4;
//         for (int i=0; i<n_nodes_elem1; i++)
//       elem1->set_node(i) = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
//
//       if (elem1->volume() < 0.0) {
//         Node * n = elem->get_node(1);
//         elem->set_node(1) = elem->get_node(2);
//         elem->set_node(2) = n;
//         n = elem->get_node(4);
//         elem->set_node(4) = elem->get_node(6);
//         elem->set_node(6) = n;
//         n = elem->get_node(9);
//         elem->set_node(9) = elem->get_node(8);
//         elem->set_node(8) = n;
//        }
//        break;
//      case INTERP_KERNEL::NORM_HEXA8:
//
//        if (elem->volume() < 0.0) {
//         Node * n = elem->get_node(1);
//         elem->set_node(1) = elem->get_node(3);
//         elem->set_node(3) = n;
//         n = elem->get_node(5);
//         elem->set_node(5) = elem->get_node(7);
//         elem->set_node(7) = n;
//
//       }
//       break;
//
//      case INTERP_KERNEL::NORM_HEXA20:
//
//
//        for (int i=0; i<n_nodes_elem; i++)
//         elem->set_node(transl[i])
//         = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
//
//
//        elem1 = new Hex8;
//        n_nodes_elem1 = 8;
//        for (int i=0; i<n_nodes_elem1; i++)
//       elem1->set_node(transl[i]) = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
//
//
//
//
//        if (elem1->volume() < 0.0) {
//
//
//         Node * n = elem->get_node(1);
//         elem->set_node(1) = elem->get_node(3);
//         elem->set_node(3) = n;
//         n = elem->get_node(8);
//         elem->set_node(8) = elem->get_node(11);
//         elem->set_node(11) = n;
//         n = elem->get_node(9);
//         elem->set_node(9) = elem->get_node(10);
//         elem->set_node(10) = n;
//         n = elem->get_node(5);
//         elem->set_node(5) = elem->get_node(7);
//         elem->set_node(7) = n;
//         n = elem->get_node(12);
//         elem->set_node(12) = elem->get_node(19);
//         elem->set_node(19) = n;
//         n = elem->get_node(13);
//         elem->set_node(13) = elem->get_node(18);
//         elem->set_node(18) = n;
//         n = elem->get_node(14);
//         elem->set_node(14) = elem->get_node(17);
//         elem->set_node(17) = n;
//
//       }
//
//       break;
//
//        case INTERP_KERNEL::NORM_HEXA27:
//
//
//        for (int i=0; i<n_nodes_elem; i++)
//         elem->set_node(transl27[i])
//         = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
//
//
//
//        elem1= new Hex8;
//        n_nodes_elem1 = 8;
//        for (int i=0; i<n_nodes_elem1; i++)
//       elem1->set_node(transl27[i]) = node_ptr(_node_id[meshConnect->getIJ(cstart+1, i)]);
//
//        if (elem1->volume() < 0.0) {
//         Node * n = elem->get_node(1);
//         elem->set_node(1) = elem->get_node(3);
//         elem->set_node(3) = n;
//         n = elem->get_node(8);
//         elem->set_node(8) = elem->get_node(11);
//         elem->set_node(11) = n;
//         n = elem->get_node(9);
//         elem->set_node(9) = elem->get_node(10);
//         elem->set_node(10) = n;
//         n = elem->get_node(5);
//         elem->set_node(5) = elem->get_node(7);
//         elem->set_node(7) = n;
//         n = elem->get_node(12);
//         elem->set_node(12) = elem->get_node(19);
//         elem->set_node(19) = n;
//         n = elem->get_node(13);
//         elem->set_node(13) = elem->get_node(18);
//         elem->set_node(18) = n;
//         n = elem->get_node(14);
//         elem->set_node(14) = elem->get_node(17);
//         elem->set_node(17) = n;
//         n = elem->get_node(22);
//         elem->set_node(22) = elem->get_node(23);
//         elem->set_node(23) = n;
//         n = elem->get_node(24);
//         elem->set_node(24) = elem->get_node(21);
//         elem->set_node(21) = n;
//
//       }
//
//       break;
//     default:
//       break;
//     }
//
//
//     elem->set_id(_elem_id[iElem]);
//
//     elem->processor_id() = iProc;
//
//     this->add_elem (elem);
//
//     if (FDEBUG) {
//       fDebug << "element (";
//       for (int i=0; i<n_nodes_elem; i++)
// 	fDebug << std::setw(6) << elem->get_node(i)->id() << " ";
//       fDebug     << "), id = " << std::setw(6) << _elem_id[iElem]
// 		 << " _proc_id " << std::setw(6) << iProc << std::endl;
//     }
//
//   }
//
//   delete [] sOut;
//   fDebug << "Elapsed     (insert nodes-elems)    " << T.elapsed() << std::endl;
//
//   _n_nodes = nodeIdOffset[nProcs];
//   _n_elem = elemIdOffset[nProcs];
//   if (FDEBUG) {
//     fDebug << "n_nodes : " << _n_nodes << " " << this->n_local_nodes() << std::endl;
//     fDebug << "n_nodes : " << this->n_nodes() << " " << this->n_local_nodes() << std::endl;
//   }
//
//   _is_serial = false;
//   _partitioner = AutoPtr<libMesh::Partitioner>(NULL);
//   T.start();
//   //  MeshCommunication().gather_neighboring_elements(*this);
//   prepare_for_use(false);
//   fDebug << "Elapsed     (prepare)    " << T.elapsed() << std::endl;
//
//   //  fDebug << boundary_info->n_boundary_ids() << std::endl;
//
//   int nBdyFaces = 0;
//   MeshBase::const_element_iterator       el
//     = active_local_elements_begin();
//   const MeshBase::const_element_iterator end_el
//     = active_local_elements_end();
//   for ( ; el != end_el; ++el) {
//     const Elem* elem = *el;
//     for (unsigned int side=0; side<elem->n_sides(); side++)
//       if (elem->neighbor(side) == NULL)
// 	{
// 	  nBdyFaces++;
// 	  if (FDEBUG) fDebug << " elem " << elem->id()
// 		 << " side : " << side << std::endl;
// 	}
//   }
//
//   delete [] nodeIdOffset;
//   delete [] elemIdOffset;
  return;
}
// =======================================================
ParallelMeshExtended::~ParallelMeshExtended() {
//   _meshes.clear();
  _node_id.clear();
  _proc_id.clear();
  _elem_id.clear();
}

// ========================================================
// This function print the med file
void ParallelMeshExtended::print_med(int Level, std::string filename) {
  
  // med->libmesh map  second order only
#if ELTYPE==27
#if DIMENSION==3
  const unsigned int nodesinv[]   = {4,7,3,0,5,6,2,1,19,15,11,12,17,14,9,13,
                                     16,18,10,8,24,25,23,20,21,22,26};
  const unsigned int nodesinvbd[] = {3,0,1,2,7,4,5,6,8};
#else
  const unsigned int nodesinv[] = {3,0,1,2,7,4,5,6,8};
  const unsigned int nodesinvbd[] = {0,2,1};
#endif
#endif
#if ELTYPE==10
#if DIMENSION==3
  const unsigned int nodesinv[]   = {2,3,1,0,9,8,5,6,7,4};
  const unsigned int nodesinvbd[] = {1,2,0,4,5,3};
#else
  const unsigned int nodesinv[]   = {1,2,0,4,5,3};
  const unsigned int nodesinvbd[] = {0,2,1};
#endif
#endif

  int       el_neigh[NDOF_FEM];          // bd element connectivity
  int       sur_toply[NDOF_FEMB];
  int       el_conn[NDOF_FEM];
  int       elb_conn[NDOF_FEMB];
  double    *xx_qnds=new double[_dim*NDOF_FEM];

  // setup  -------------------------------------------
//   Level=_NoLevels-1;
  int      n_nodes    =_NoNodes[Level];
  int      n_elements = _NoElements[0][Level];
  int       nodes_el   = _type_FEM[0];
  const int el_sides=  _GeomEl._n_sides[0];
  const int pt_sides_bd=NDOF_FEMB;
//
//   // init  -------------------------------------------
  std::cout << "n_nodes " << n_nodes << "n_elements "
            << n_elements << "nodes_el " << nodes_el << std::endl;
 int icount =0;  int ielemcount =0;
  int ibccount=0; int bd_n_elements =0;

  _node_id.resize(n_nodes);
  _proc_id.resize(n_nodes);
  _elem_id.resize(n_elements);
//
//   // coords ------------------------------------------
  double * coord; coord = new double[n_nodes*_dim];
  for(int i=0; i<n_nodes; i++) {
    for(int idim=0; idim<_dim; idim++) {
      coord[i*_dim+idim]=_xyz[i+idim*n_nodes];
    }
    _node_id[i]=i+1; // map my-mesh->MEDcoupling
  }
//
//
//
  int * conn_bd; conn_bd=new int [n_elements*el_sides*pt_sides_bd];
  int * conn; conn=new int [n_elements*nodes_el];
  std::vector<int> elem_bd_id2;
//   elem_id2[0]=0;
//   std::vector< int > asf; asf.resize(n_elements);
  for(int  iproc = 0; iproc <_n_subdom; iproc++) {
    const int nel_e = _off_el[0][Level+_NoLevels*iproc+1]; // start element
    const int nel_b = _off_el[0][Level+_NoLevels*iproc];   // stop element
    for(int el=0; el < (nel_e - nel_b); el++) {
//     for (int el =0 ;
//          el <_off_el[0][iproc*_NoLevels+Level+1]-
//          _off_el[0][iproc*_NoLevels+Level]; el++) {
      get_el_nod_conn(0,Level,el,el_conn,xx_qnds,iproc);
      get_el_neighbor(el_sides,0,Level,el,el_neigh,iproc);

      for(int  i = 0; i < nodes_el; i++) {
        conn[icount] = el_conn[nodesinv[i]];//el_map[0][el*nodes_el+i];
         _proc_id[conn[icount]]= iproc;// map my-mesh->iproc
         icount++;
      }
      // element boundary
      for(int iside=0; iside< el_sides; iside++)  {
        if(el_neigh[iside] == -1){   // domain boundary
          int min=1000000;
          for(int idof=0; idof<pt_sides_bd; idof++) {                          //idof -> boundary node
            sur_toply[idof]=_GeomEl._surf_top[nodesinvbd[idof]+pt_sides_bd*iside]; //use map to find global node
            int idofb=sur_toply[idof];   //idofb -> element node
            conn_bd[ibccount]=el_conn[idofb];//_el_map[0][el*nodes_el+idofb];
            ibccount++;
//
            if(min>fabs(_bc_id[el_conn[idofb]]))  min=fabs(_bc_id[el_conn[idofb]]);
          }
          elem_bd_id2.push_back(min);
          bd_n_elements++;  // domain boundary found
        }
      }
      _elem_id[el]=ielemcount+1; // map my-mesh->MEDcoupling

//       asf[ielemcount]=el+nel_b;

      ielemcount++;
    }  // end for-el
  } // end for-iproc
  

  // MED mesh *************************************************
//
  // MEDCouplingUMesh mesh connectivity (volume)
  ParaMEDMEM::MEDCouplingUMesh *mesh1=ParaMEDMEM::MEDCouplingUMesh::New("Mesh_1",_dim);
  mesh1->allocateCells(n_elements);
  for(int  i = 0; i < n_elements; i++)
    mesh1->insertNextCell(MED_EL_TYPE,nodes_el,conn+i*nodes_el);
  mesh1->finishInsertingCells();

// MEDCouplingUMesh Mesh connectivity (boundary)
  ParaMEDMEM::MEDCouplingUMesh *mesh2=ParaMEDMEM::MEDCouplingUMesh::New("Mesh_1",DIMENSION-1);
  mesh2->allocateCells(bd_n_elements);
  for(int  i = 0; i < bd_n_elements; i++)
    mesh2->insertNextCell(MED_EL_BDTYPE,pt_sides_bd,conn_bd+i*pt_sides_bd);
  mesh2->finishInsertingCells();

  // coord (same node set for both meshes)
  ParaMEDMEM::DataArrayDouble *coordarr=ParaMEDMEM::DataArrayDouble::New();
  coordarr->alloc(n_nodes,_dim);
  std::copy(coord,coord+n_nodes*_dim,coordarr->getPointer());
  mesh1->setCoords(coordarr);
  mesh2->setCoords(coordarr);
// //
  // Setting MEDCouplingUMesh into MEDFileUMesh
  ParaMEDMEM::MEDFileUMesh *mm=ParaMEDMEM::MEDFileUMesh::New();
  mm->setName("Mesh_1");//name needed to be non empty
  mm->setDescription("Description Mesh_1");
  mm->setCoords(mesh1->getCoords());
  mm->setMeshAtLevel(0,mesh1,false);
  mm->setMeshAtLevel(-1,mesh2,false);

  //Volume Groups

  std::map <int,std::vector<int> > vol_group_elements;
  for(int i=0; i<n_elements; i++) {
    vol_group_elements[_mat_id[i]].push_back(i);
  }

  int n_vol_group=vol_group_elements.size();
  std::vector<const ParaMEDMEM::DataArrayInt *> gr_vol(n_vol_group);
  ParaMEDMEM::DataArrayInt ** g_vol=new ParaMEDMEM::DataArrayInt *[n_vol_group];

  int js=0;
  // defining the vol group data to store
  std::map<int,std::vector<int> >::iterator it_vol;
  for(it_vol=vol_group_elements.begin(); it_vol!=vol_group_elements.end(); ++it_vol) {
    int igroup=it_vol->first;  int is=it_vol->second.size();
    g_vol[js]=ParaMEDMEM::DataArrayInt::New();
    g_vol[js]->alloc(is,1);
    std::ostringstream name_p; name_p<< igroup;
    g_vol[js]->setName(name_p.str().c_str());
    int *val1=new int[is];
    for(int iv=0; iv<is; iv++)  val1[iv]=it_vol->second[iv];
    std::copy(val1,val1+is,g_vol[js]->getPointer());
    delete []val1;
    gr_vol[js]=g_vol[js] ;
    js++;
  }
  // inserting  the volume groups into the med-file
  mm->setGroupsAtLevel(0,gr_vol,false);


  // Boundary group ******************************************
  // Finding bd_group and  bd_group_elements map
  std::map <int,std::vector<int> > bd_group_elements;
//   std::map <int, int> bd_group;
  for(int i=0; i<bd_n_elements ; i++) {
//     bd_group [elem_bd_id2[i]]++;
    bd_group_elements[elem_bd_id2[i]].push_back(i);
  }
  int n_bd_group=bd_group_elements.size();
  // group vector
  std::vector<const ParaMEDMEM::DataArrayInt *> gr_bd(n_bd_group);
  ParaMEDMEM::DataArrayInt ** g_bd=new ParaMEDMEM::DataArrayInt *[n_bd_group];

  js=0;
  // defining the  group data to store
  std::map<int,std::vector<int> >::iterator it;
  for(it=bd_group_elements.begin(); it!=bd_group_elements.end(); ++it) {
    int igroup=it->first;  int is=it->second.size();
//     std::cout << igroup << "\n";
    g_bd[js]=ParaMEDMEM::DataArrayInt::New();
    g_bd[js]->alloc(is,1);
    std::ostringstream name_p; name_p<< igroup;
    g_bd[js]->setName(name_p.str().c_str());
    int *valb1=new int[is]; /*int icount=0;*/
    for(int iv=0; iv<is; iv++)  valb1[iv]=it->second[iv];
    std::copy(valb1,valb1+is,g_bd[js]->getPointer());
    delete []valb1;
    gr_bd[js]=g_bd[js] ;
    js++;
  }
  // insert the boundary groups into the med-file
  mm->setGroupsAtLevel(-1,gr_bd,false);

  // Printing Group and Family
  std::cout << "\n \n =================================== ";
  std::cout << "\n Group Names -> Families : \n";
  std::map<std::string, std::vector<std::string> > a=(mm->getGroupInfo());
  std::map<std::string, std::vector<std::string> >::iterator ita;
  for(ita=a.begin();  ita!=a.end(); ++ita) {
    std::string igroup=ita->first;  int is=ita->second.size();
    std::cout << "\n "<< igroup << "-> ";
    for(int i=0; i<is; i++) {
      std::string a2 =ita->second[i];
      std::cout << a2 << "  ";
    }
  }
  std::cout << "\n \n =================================== ";
  std::cout << " \n Families -> Groups id: \n";
  std::map<std::string,int> fama=mm->getFamilyInfo();
  std::map<std::string, int >::iterator itfa;
  for(itfa=fama.begin();  itfa!=fama.end(); ++itfa) {
    std::string igroup=itfa->first;
    std::cout << "\n "<< igroup << "-> ";
    int a2 =itfa->second;
    std::cout << a2 << "  ";
  }
  //mm->write(filename.c_str(),2);

  // clean
  delete []coord; fama.clear(); a.clear();
   coordarr->decrRef(); mesh1->decrRef(); mesh2->decrRef();
   mm->decrRef();
   delete []conn_bd;   delete []conn;
   elem_bd_id2.clear();
   vol_group_elements.clear();bd_group_elements.clear();
return;
}

void ParallelMeshExtended::read_bc_id(int Level) {

  std::string    input_dir = _mgutils._inout_dir;
  std::string    basemesh = _mgutils.get_file("BASEMESH");

  // Open an existing file. ---------------
  std::ostringstream meshname;  meshname << input_dir << basemesh << ".h5";
  std::cout << " Reading bc_id from= " <<  meshname.str() <<  std::endl;
  hid_t  file_id = H5Fopen(meshname.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

  // Reading  bc_id ---------------------------------------

  // Getting dataset
  std::ostringstream Name; Name << "NODES/COORD/BC";
  hid_t dtset = H5Dopen(file_id,Name.str().c_str()
#if HDF5_VERSIONM!=188
                        ,H5P_DEFAULT
#endif
                       );
  hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  hsize_t dims[2];
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status ==0) std::cerr << "SalomeIO::read dims not found";
//   _n_bc_id=dims[0];
  _bc_id.resize(dims[0]);

  // reading
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&_bc_id[0]);
  
  
  //   int ielemn=0;
  for (int iproc=0; iproc<_n_subdom; iproc++) {
    int nel_b=_off_el[0][_NoLevels*iproc+Level];
    int nel_e=_off_el[0][Level+_NoLevels*iproc+1];
    for (int iel=nel_b; iel <nel_e; iel++) {
      
     
      for (int  inode=NDOF_P; inode<NDOF_FEM; inode++)    {
        //get the global node number
        int el_conn = _el_map[0][iel*NDOF_FEM+inode];
        // get the element coordinades
        _bc_id[el_conn] *=-1;
      }
  } // ---- end iel -------------------------
} // 2bB end interpolation over the fine mesh -
  
  
  
  
  
  return;
}

void ParallelMeshExtended::read_mat_id(int Level) {

  std::string    input_dir = _mgutils._inout_dir;
  std::string    basemesh = _mgutils.get_file("BASEMESH");

  // Open an existing file. ---------------
  std::ostringstream meshname;  meshname << input_dir << basemesh << ".h5";
  std::cout << " Reading mat_id from= " <<  meshname.str() <<  std::endl;
  hid_t  file_id = H5Fopen(meshname.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

  // Reading  bc_id ---------------------------------------

  // Getting dataset
  std::ostringstream Name; Name << "ELEMS/SUB/MAT"<<Level;
  hid_t dtset = H5Dopen(file_id,Name.str().c_str()
#if HDF5_VERSIONM!=188
                        ,H5P_DEFAULT
#endif
                       );
  hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  hsize_t dims[2];
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status <0) {std::cerr << " ******** SalomeIO::read dims not found"; abort();}
   _mat_id.resize(dims[0]);

  // reading
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&_mat_id[0]);

  return;
}
