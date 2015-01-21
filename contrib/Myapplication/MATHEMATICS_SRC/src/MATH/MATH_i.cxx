
// this cxx file was generated by yacsgen
#include "MATH_i.hxx"
#include "MATH.hxx"
using namespace std;
#include <string>
#include <vector>
#include "SenderFactory.hxx"
#include "MultiCommException.hxx"
#include "ReceiverFactory.hxx"
#include "SALOME_Matrix_i.hxx"
#include "MatrixClient.hxx"
#include "Utils_CorbaException.hxx"

//DEFS






//ENDDEF

//=============================================================================
/*!
 *  standard constructor
 */
//=============================================================================
MATH_i::MATH_i(CORBA::ORB_ptr orb,
	PortableServer::POA_ptr poa,
	PortableServer::ObjectId * contId, 
	const char *instanceName, 
	const char *interfaceName) :
  Engines_Component_i(orb, poa, contId, instanceName, interfaceName),cppCompo_(new MATH)
{
  MESSAGE("activate object");
  _thisObj = this ;
  _id = _poa->activate_object(_thisObj);
}

MATH_i::~MATH_i()
{
}


CORBA::Double MATH_i::sum(const MATHEMATICS_ORB::dblevec& tab) throw (SALOME::SALOME_Exception)
{
    beginService("MATH_i::sum");
    BEGIN_OF("MATH_i::sum");
    try
    {
//	Arguments processing
	long _tab_size=tab.length();
	const double *_tab_value = &tab[0];
	std::vector<double> _tab(_tab_value,_tab_value+_tab_size);
//	Call cpp component
	double _rtn_cpp = cppCompo_->sum( _tab);
//	Post-processing & return
	CORBA::Double _rtn_ior(_rtn_cpp);
	endService("MATH_i::sum");
	END_OF("MATH_i::sum");
	return _rtn_ior;
    }
    catch (std::exception& ex)
    {
        THROW_SALOME_CORBA_EXCEPTION( ex.what(), SALOME::INTERNAL_ERROR );
    }
}


CORBA::Double MATH_i::squareroot(CORBA::Double x) throw (SALOME::SALOME_Exception)
{
    beginService("MATH_i::squareroot");
    BEGIN_OF("MATH_i::squareroot");
    try
    {
//	Arguments processing
	double _x(x);
//	Call cpp component
	double _rtn_cpp = cppCompo_->squareroot( _x);
//	Post-processing & return
	CORBA::Double _rtn_ior(_rtn_cpp);
	endService("MATH_i::squareroot");
	END_OF("MATH_i::squareroot");
	return _rtn_ior;
    }
    catch (std::exception& ex)
    {
        THROW_SALOME_CORBA_EXCEPTION( ex.what(), SALOME::INTERNAL_ERROR );
    }
}


CORBA::Double MATH_i::plus(CORBA::Double x,CORBA::Double y) throw (SALOME::SALOME_Exception)
{
    beginService("MATH_i::plus");
    BEGIN_OF("MATH_i::plus");
    try
    {
//	Arguments processing
	double _x(x);
	double _y(y);
//	Call cpp component
	double _rtn_cpp = cppCompo_->plus( _x, _y);
//	Post-processing & return
	CORBA::Double _rtn_ior(_rtn_cpp);
	endService("MATH_i::plus");
	END_OF("MATH_i::plus");
	return _rtn_ior;
    }
    catch (std::exception& ex)
    {
        THROW_SALOME_CORBA_EXCEPTION( ex.what(), SALOME::INTERNAL_ERROR );
    }
}


CORBA::Double MATH_i::minus(CORBA::Double x,CORBA::Double y) throw (SALOME::SALOME_Exception)
{
    beginService("MATH_i::minus");
    BEGIN_OF("MATH_i::minus");
    try
    {
//	Arguments processing
	double _x(x);
	double _y(y);
//	Call cpp component
	double _rtn_cpp = cppCompo_->minus( _x, _y);
//	Post-processing & return
	CORBA::Double _rtn_ior(_rtn_cpp);
	endService("MATH_i::minus");
	END_OF("MATH_i::minus");
	return _rtn_ior;
    }
    catch (std::exception& ex)
    {
        THROW_SALOME_CORBA_EXCEPTION( ex.what(), SALOME::INTERNAL_ERROR );
    }
}


CORBA::Double MATH_i::times(CORBA::Double x,CORBA::Double y) throw (SALOME::SALOME_Exception)
{
    beginService("MATH_i::times");
    BEGIN_OF("MATH_i::times");
    try
    {
//	Arguments processing
	double _x(x);
	double _y(y);
//	Call cpp component
	double _rtn_cpp = cppCompo_->times( _x, _y);
//	Post-processing & return
	CORBA::Double _rtn_ior(_rtn_cpp);
	endService("MATH_i::times");
	END_OF("MATH_i::times");
	return _rtn_ior;
    }
    catch (std::exception& ex)
    {
        THROW_SALOME_CORBA_EXCEPTION( ex.what(), SALOME::INTERNAL_ERROR );
    }
}


CORBA::Double MATH_i::sinx(CORBA::Double x) throw (SALOME::SALOME_Exception)
{
    beginService("MATH_i::sinx");
    BEGIN_OF("MATH_i::sinx");
    try
    {
//	Arguments processing
	double _x(x);
//	Call cpp component
	double _rtn_cpp = cppCompo_->sinx( _x);
//	Post-processing & return
	CORBA::Double _rtn_ior(_rtn_cpp);
	endService("MATH_i::sinx");
	END_OF("MATH_i::sinx");
	return _rtn_ior;
    }
    catch (std::exception& ex)
    {
        THROW_SALOME_CORBA_EXCEPTION( ex.what(), SALOME::INTERNAL_ERROR );
    }
}


char* MATH_i::ComponentDataType()
{
      return CORBA::string_dup("MATH");
}

Engines::EngineComponent_ptr MATH_i::GetComponentInstance()
{
      return MATH_Gen::_this();
}

extern "C"
{
  PortableServer::ObjectId * MATHEngine_factory(
			       CORBA::ORB_ptr orb,
			       PortableServer::POA_ptr poa, 
			       PortableServer::ObjectId * contId,
			       const char *instanceName, 
		       	       const char *interfaceName)
  {
    MESSAGE("PortableServer::ObjectId * MATHEngine_factory()");
    SCRUTE(interfaceName);
    MATH_i * myMATH 
      = new MATH_i(orb, poa, contId, instanceName, interfaceName);
    return myMATH->getId() ;
  }
}
