
// class include
#include "conf.hxx"
#include "FEMUS.hxx"

// MED includes
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "InterpKernelExprParser.hxx"

// LibMeshCpp includes
#include "ParallelMeshExtended.hxx"

#include "_LibMeshProblem.hxx"
#include "Debug.hxx"

// Libmesh test case includes
#include "LaplaceProblem.hxx"
#include "NS_Problem.hxx"
#include "NS_Problem3D.hxx"
#include "LinearElasticProblem.hxx"
#include "Con_DiffProblem.hxx"
#include "MonodimProblem.hxx"

// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// ============================================================================
// Basic constructor
FEMUS::FEMUS() :
  _comm(MPI_COMM_WORLD),
  _problem(NULL),
  _mesh(NULL) { //=============================================================
  int flag;
  MPI_Initialized(&flag);

  if (flag)  _local_MPI_Init = false;
  else {
    _local_MPI_Init = true;
    int argc = 0; char ** argv = NULL;
    MPI_Init(&argc, &argv);
  }
}
// ============================================================================
// This function is a constructor with  communicator
FEMUS::FEMUS(
  MPI_Comm comm
) :
  _comm(comm),
  _problem(NULL),
  _mesh(NULL) { // ============================================================
  int flag;
  MPI_Initialized(&flag);

  if (flag)  _local_MPI_Init = false;
  else {
    _local_MPI_Init = true;
    int argc = 0; char ** argv = NULL;
    MPI_Init(&argc, &argv);
  }
}
// ============================================================================
// This function is the problem destructor
void FEMUS::terminate(
) {// =========================================================================
  if (_problem) _problem->terminate();
}

// ============================================================================
// This function is the destructor
FEMUS::~FEMUS() {
  // ==========================================================================
  if (_problem) delete _problem;     // problem
  if (_mesh) _mesh->decrRef();       // med-mesh
  for (int i=0; i<_bc.size(); i++)  if (_bc[i].support) _bc[i].support->decrRef();
  if (_local_MPI_Init)  MPI_Finalize();
}
// ****************************************************************************
// ****************    end Constructor Destructor *****************************

// ****************************************************************************
// ****************    Set    *************************************************

// =============================================================================
// This function sets the type of problem
// Poisson
// Elasticity
// convdiff     monodim
// navierstokes navierstokes3D
void FEMUS::setType(
  const std::string & pbName
) {// ==========================================================================
  std::string sPbName = pbName;

  if (sPbName == "Poisson")                        // Laplace 1D 2D 3D
    _problem = new LaplaceProblem(_comm);
  else if (sPbName == "Elasticity")                // Elaticity 
    _problem = new LinearElasticProblem(_comm);
  else if (sPbName == "convdiff")                  // Heat equation
    _problem = new Con_DiffProblem(_comm);
  else if (sPbName == "navierstokes")              // Navier-Stokes 2D
    _problem = new NS_Problem(_comm);
  else if (sPbName == "navierstokes3D")            // Navier-Stokes 3D
    _problem = new NS_Problem3D(_comm);
  else if (sPbName == "monodim")                   // monodimensional
    _problem = new Monodim_Problem(_comm);
  else
    _problem = NULL;

  return;
}



// // from French speaker
// std::string replaceEnv(const std::string & s) {
//   std::string t = s;
//   size_t p0 = t.find_first_of('$'), p1 = std::string::npos;
//   if (p0 != std::string::npos) {
//     p1 = t.find_first_not_of
//          ("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789",p0+1);
//     if (p1 != std::string::npos) p1--;
// 
//     std::string env = t.substr(p0+1, p1-p0);
//     char * pEnv = getenv(env.c_str());
//     if (pEnv != NULL) {
//       t = t.replace(p0, p1+1, pEnv);
//     }
//   }
//   return t;
// }


// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void FEMUS::setMesh(
  const ParaMEDMEM::MEDCouplingUMesh * m
) {
  if (_problem) _problem->setMesh(m);
}

// =============================================================================
// This function reads and sets the med-mesh 
// (also the libmesh calling the other setMesh function)
void FEMUS::setMesh(const std::string & dataFile) {
  
  int iProc; MPI_Comm_rank(_comm, &iProc);

//   std::string t = replaceEnv(dataFile);
  // Reading MED mesh ----------------------------------------------------
  std::string localFile; int l = dataFile.size();
  if (dataFile.substr(l-4) == ".med")    localFile = dataFile;
  else {
    std::ostringstream s1;
    s1 << dataFile << iProc+1 << ".med";
    localFile = s1.str();
  }
  std::cerr << "MED file for process " << iProc << ": "<< localFile << std::endl;

  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
  if (meshNames.size() < 1) {
    std::string message = "data file '";
    message += localFile;
    message += "' doesn't contain a mesh (in MED format)";
    throw std::runtime_error(message);
  }

  std::string localMeshName = meshNames[0];
  _mesh = MEDLoader::ReadUMeshFromFile(localFile.c_str(), localMeshName.c_str(), 0);
  
  if (_mesh == NULL) {
    std::string message = "unable to read '";
    message += localMeshName;
    message += "' from file '";
    message += dataFile;
    throw std::runtime_error(message);
  }
  
  // set FEMUS mesh ----------------------------------------
  setMesh(_mesh);

  // Reading boundary conditions ------------------------------
  std::vector<std::string> GroupNames = 
  MEDLoader::GetMeshGroupsNames(localFile.c_str(), localMeshName.c_str());

  int nGroups = GroupNames.size();
  std::vector<std::string> vG(1);
  _bc.resize(nGroups);
  for (int i=0; i<nGroups; i++) {

    vG[0] = GroupNames[i];
    _bc[i].name = GroupNames[i];
    _bc[i].support = MEDLoader::ReadUMeshFromGroups
                     (localFile.c_str(), meshNames[0].c_str(), -1, vG);
    if (_bc[i].support) {
      _bc[i].support->zipCoords();
      _bc[i].id = defineBoundary(_bc[i].support);
    } else
      _bc[i].id = -1;
  }
  
  return;
}

// *******************************************************************
// ********************* set and get  boundary conditions ************

// ===================================================================
std::vector<std::string> FEMUS::getBoundaryNames(
) { //================================================================
  std::vector<std::string> names(_bc.size());
  for (int i = 0; i<_bc.size(); i++)
    names[i] = _bc[i].name;
  return names;
}

std::string FEMUS::getBoundaryName(
  int i) {
  if (i>=0 && i<_bc.size())
    return _bc[i].name;
  return "";
}

void FEMUS::setBCType
(const std::string & nameBC, const std::string & typeBC) {
  int id = searchBoundaryCondition(nameBC);
  if (_problem) {
    _problem->setBoundaryConditionType
    (id, typeBC.c_str());
    _problem->setAnalyticBoundaryValues(id, "0.0");
  }
}

void FEMUS::setAnalyticBoundaryValues(const std::string & bcName,
                                        const std::string & bcType,
                                        const std::string & bcExpression) {
  if (FDEBUG) fDebugPos << bcName << " " << bcType << std::endl;
  setBCType(bcName, bcType);
  if (bcType != "Neumann") setAnalyticBCValues(bcName, bcExpression);
}


void FEMUS::setFieldBoundaryValues
(const std::string & bcName,
 const std::string & bcType,
 const ParaMEDMEM::MEDCouplingFieldDouble * bcField) {
  if (FDEBUG) fDebugPos << bcName << " " << bcType << std::endl;
  setBCType(bcName, bcType);
  setFieldBCValues(bcName, bcField);
}

void FEMUS::setAnalyticBCValues
(const std::string & nameBC, const std::string & f) {
  int id = searchBoundaryCondition(nameBC);
  if (_problem) _problem->setAnalyticBoundaryValues(id, f.c_str());
}

void FEMUS::setFieldBCValues
(const std::string & nameBC, const ParaMEDMEM::MEDCouplingFieldDouble *f) {
  int id = searchBoundaryCondition(nameBC);
  if (_problem) _problem->setFieldBoundaryValues(id, f);
}

int FEMUS::defineBoundary(const ParaMEDMEM::MEDCouplingUMesh * mesh_bc) {
  return (_problem) ? _problem->defineBoundary(mesh_bc) : -1;
}

ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary
(const std::string & nameBC, const std::string & vName) const {
  int id = searchBoundaryCondition(nameBC);
  return _problem ? _problem->getValuesOnBoundary(id, vName.c_str()) : NULL;
}

ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary
(const std::string & nameBC,  std::vector<char *> vName) const {
  int id = searchBoundaryCondition(nameBC);
  return _problem ? _problem->getValuesOnBoundary(id, vName) : NULL;
}

ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary_nodes
(const std::string & nameBC,  std::vector<char *> vName) const {
  int id = searchBoundaryCondition(nameBC);
  return _problem ? _problem->getValuesOnBoundary_nodes(id, vName) : NULL;
}

int FEMUS::searchBoundaryCondition(const std::string & nameBC) const {
  for (std::vector<sBC>::const_iterator i = _bc.begin(); i != _bc.end(); i++)
    if (i->name == nameBC) return i->id;

  std::string message = "non existent boundary condition named ";
  message += nameBC;
  message += "(choose between:";
  for (std::vector<sBC>::const_iterator i = _bc.begin(); i != _bc.end(); i++) {
    message += " '";
    message += i->name;
    message += "'";
  }
  message += ")";
  throw std::runtime_error(message);
}



// *******************************************************************
// **************** source *******************************************

// ===================================================================
// This function sets an analytic source
void FEMUS::setAnalyticSource(
  const std::string & f                       // analytic source
) { // ================================================================
  if (_problem) _problem->setAnalyticSource(f.c_str());
}

// ===================================================================
// This function sets a numerical source (coupling)
void FEMUS::setSource(
  const ParaMEDMEM::MEDCouplingFieldDouble *f  // numerical sources
) { // ================================================================
  if (_problem) _problem->setSource(f);
}

// ===================================================================
// This function sets a numerical source (coupling)
void FEMUS::setAverageSources(
  int n,              //  number of average sources <-
  double val[]        //  average source            <-
) {
  if (_problem) _problem->setAverageSources(n,val);
}

// *******************************************************************
// *******************************************************************
// **************** Solve  *******************************************
//====================================================================
// This function solves the problem
void FEMUS::solve() {
  if (_problem) _problem->solve();
}

// *******************************************************************
// **************** Field solution************************************

//====================================================================
//====================================================================
// This function gets the solution in MED format
ParaMEDMEM::MEDCouplingFieldDouble *FEMUS::getOutputField(
  const std::string & vName         // field name <-    
) const { //==========================================================
  return (_problem) ? _problem->getOutputField(vName.c_str()) : NULL;
}

//====================================================================
// This function prints the solution in MED format
void FEMUS::saveOutputField( 
  const std::string & vName,                // field name <-
  const std::string & prefix                // file name  <-
) const { // ==========================================================
  ParaMEDMEM::MEDCouplingFieldDouble * f = getOutputField(vName);
//   std::string t = replaceEnv(prefix);

  int iProc;  MPI_Comm_rank(_comm, &iProc);
  std::ostringstream s;
  s << prefix << "_" << vName << "_" << iProc+1 << ".med";
  MEDLoader::WriteField(s.str().c_str(), f, true);
  MPI_Barrier(_comm);
}




//====================================================================
// This function print the solution in general format
ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getInputFieldTemplate(const std::string &name) {
//    if (name == "")
//       return getOutputField("???");
//    else
//       return getValuesOnBoundary(name, "???");
}



// *******************************************************************
// **************** Debug *******************************************

// ===================================================================
// This function sets the file for logging
void FEMUS::setDebugFile(
  const std::string & s
) { // ===============================================================
  if (FDEBUG) fDebug.open(s.c_str(), _comm);
}

// ===================================================================
// This function ends the  logging
void FEMUS::endDebug(
) {// ===============================================================
  fDebug.close();
}
