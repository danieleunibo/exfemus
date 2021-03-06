// This file is generated by omniidl (C++ backend) - omniORB_4_1. Do not edit.

#include "MATHEMATICS.hh"

OMNI_USING_NAMESPACE(omni)

static const char* _0RL_dyn_library_version = omniORB_4_1_dyn;

static ::CORBA::TypeCode::_Tracker _0RL_tcTrack(__FILE__);

static CORBA::TypeCode_ptr _0RL_tc_MATHEMATICS__ORB_mstringvec = CORBA::TypeCode::PR_alias_tc("IDL:MATHEMATICS_ORB/stringvec:1.0", "stringvec", CORBA::TypeCode::PR_sequence_tc(0, CORBA::TypeCode::PR_string_tc(0, &_0RL_tcTrack), &_0RL_tcTrack), &_0RL_tcTrack);


#if defined(HAS_Cplusplus_Namespace) && defined(_MSC_VER)
// MSVC++ does not give the constant external linkage otherwise.
namespace MATHEMATICS_ORB { 
  const ::CORBA::TypeCode_ptr _tc_stringvec = _0RL_tc_MATHEMATICS__ORB_mstringvec;
} 
#else
const ::CORBA::TypeCode_ptr MATHEMATICS_ORB::_tc_stringvec = _0RL_tc_MATHEMATICS__ORB_mstringvec;
#endif

static CORBA::TypeCode_ptr _0RL_tc_MATHEMATICS__ORB_mdblevec = CORBA::TypeCode::PR_alias_tc("IDL:MATHEMATICS_ORB/dblevec:1.0", "dblevec", CORBA::TypeCode::PR_sequence_tc(0, CORBA::TypeCode::PR_double_tc(), &_0RL_tcTrack), &_0RL_tcTrack);


#if defined(HAS_Cplusplus_Namespace) && defined(_MSC_VER)
// MSVC++ does not give the constant external linkage otherwise.
namespace MATHEMATICS_ORB { 
  const ::CORBA::TypeCode_ptr _tc_dblevec = _0RL_tc_MATHEMATICS__ORB_mdblevec;
} 
#else
const ::CORBA::TypeCode_ptr MATHEMATICS_ORB::_tc_dblevec = _0RL_tc_MATHEMATICS__ORB_mdblevec;
#endif

static CORBA::TypeCode_ptr _0RL_tc_MATHEMATICS__ORB_mintvec = CORBA::TypeCode::PR_alias_tc("IDL:MATHEMATICS_ORB/intvec:1.0", "intvec", CORBA::TypeCode::PR_sequence_tc(0, CORBA::TypeCode::PR_long_tc(), &_0RL_tcTrack), &_0RL_tcTrack);


#if defined(HAS_Cplusplus_Namespace) && defined(_MSC_VER)
// MSVC++ does not give the constant external linkage otherwise.
namespace MATHEMATICS_ORB { 
  const ::CORBA::TypeCode_ptr _tc_intvec = _0RL_tc_MATHEMATICS__ORB_mintvec;
} 
#else
const ::CORBA::TypeCode_ptr MATHEMATICS_ORB::_tc_intvec = _0RL_tc_MATHEMATICS__ORB_mintvec;
#endif

static CORBA::PR_structMember _0RL_structmember_Engines_mdataref[] = {
  {"ref", CORBA::TypeCode::PR_string_tc(0, &_0RL_tcTrack)}
};

#ifdef _0RL_tc_Engines_mdataref
#  undef _0RL_tc_Engines_mdataref
#endif
static CORBA::TypeCode_ptr _0RL_tc_Engines_mdataref = CORBA::TypeCode::PR_struct_tc("IDL:Engines/dataref:1.0", "dataref", _0RL_structmember_Engines_mdataref, 1, &_0RL_tcTrack);




static CORBA::TypeCode_ptr _0RL_tc_MATHEMATICS__ORB_mdataref = CORBA::TypeCode::PR_alias_tc("IDL:MATHEMATICS_ORB/dataref:1.0", "dataref", _0RL_tc_Engines_mdataref, &_0RL_tcTrack);


#if defined(HAS_Cplusplus_Namespace) && defined(_MSC_VER)
// MSVC++ does not give the constant external linkage otherwise.
namespace MATHEMATICS_ORB { 
  const ::CORBA::TypeCode_ptr _tc_dataref = _0RL_tc_MATHEMATICS__ORB_mdataref;
} 
#else
const ::CORBA::TypeCode_ptr MATHEMATICS_ORB::_tc_dataref = _0RL_tc_MATHEMATICS__ORB_mdataref;
#endif

#if defined(HAS_Cplusplus_Namespace) && defined(_MSC_VER)
// MSVC++ does not give the constant external linkage otherwise.
namespace MATHEMATICS_ORB { 
  const ::CORBA::TypeCode_ptr _tc_MATH_Gen = CORBA::TypeCode::PR_interface_tc("IDL:MATHEMATICS_ORB/MATH_Gen:1.0", "MATH_Gen", &_0RL_tcTrack);
} 
#else
const ::CORBA::TypeCode_ptr MATHEMATICS_ORB::_tc_MATH_Gen = CORBA::TypeCode::PR_interface_tc("IDL:MATHEMATICS_ORB/MATH_Gen:1.0", "MATH_Gen", &_0RL_tcTrack);
#endif

static void _0RL_MATHEMATICS__ORB_mstringvec_marshal_fn(cdrStream& _s, void* _v)
{
  MATHEMATICS_ORB::stringvec* _p = (MATHEMATICS_ORB::stringvec*)_v;
  *_p >>= _s;
}
static void _0RL_MATHEMATICS__ORB_mstringvec_unmarshal_fn(cdrStream& _s, void*& _v)
{
  MATHEMATICS_ORB::stringvec* _p = new MATHEMATICS_ORB::stringvec;
  *_p <<= _s;
  _v = _p;
}
static void _0RL_MATHEMATICS__ORB_mstringvec_destructor_fn(void* _v)
{
  MATHEMATICS_ORB::stringvec* _p = (MATHEMATICS_ORB::stringvec*)_v;
  delete _p;
}

void operator<<=(::CORBA::Any& _a, const MATHEMATICS_ORB::stringvec& _s)
{
  MATHEMATICS_ORB::stringvec* _p = new MATHEMATICS_ORB::stringvec(_s);
  _a.PR_insert(_0RL_tc_MATHEMATICS__ORB_mstringvec,
               _0RL_MATHEMATICS__ORB_mstringvec_marshal_fn,
               _0RL_MATHEMATICS__ORB_mstringvec_destructor_fn,
               _p);
}
void operator<<=(::CORBA::Any& _a, MATHEMATICS_ORB::stringvec* _sp)
{
  _a.PR_insert(_0RL_tc_MATHEMATICS__ORB_mstringvec,
               _0RL_MATHEMATICS__ORB_mstringvec_marshal_fn,
               _0RL_MATHEMATICS__ORB_mstringvec_destructor_fn,
               _sp);
}

::CORBA::Boolean operator>>=(const ::CORBA::Any& _a, MATHEMATICS_ORB::stringvec*& _sp)
{
  return _a >>= (const MATHEMATICS_ORB::stringvec*&) _sp;
}
::CORBA::Boolean operator>>=(const ::CORBA::Any& _a, const MATHEMATICS_ORB::stringvec*& _sp)
{
  void* _v;
  if (_a.PR_extract(_0RL_tc_MATHEMATICS__ORB_mstringvec,
                    _0RL_MATHEMATICS__ORB_mstringvec_unmarshal_fn,
                    _0RL_MATHEMATICS__ORB_mstringvec_marshal_fn,
                    _0RL_MATHEMATICS__ORB_mstringvec_destructor_fn,
                    _v)) {
    _sp = (const MATHEMATICS_ORB::stringvec*)_v;
    return 1;
  }
  return 0;
}

static void _0RL_MATHEMATICS__ORB_mdblevec_marshal_fn(cdrStream& _s, void* _v)
{
  MATHEMATICS_ORB::dblevec* _p = (MATHEMATICS_ORB::dblevec*)_v;
  *_p >>= _s;
}
static void _0RL_MATHEMATICS__ORB_mdblevec_unmarshal_fn(cdrStream& _s, void*& _v)
{
  MATHEMATICS_ORB::dblevec* _p = new MATHEMATICS_ORB::dblevec;
  *_p <<= _s;
  _v = _p;
}
static void _0RL_MATHEMATICS__ORB_mdblevec_destructor_fn(void* _v)
{
  MATHEMATICS_ORB::dblevec* _p = (MATHEMATICS_ORB::dblevec*)_v;
  delete _p;
}

void operator<<=(::CORBA::Any& _a, const MATHEMATICS_ORB::dblevec& _s)
{
  MATHEMATICS_ORB::dblevec* _p = new MATHEMATICS_ORB::dblevec(_s);
  _a.PR_insert(_0RL_tc_MATHEMATICS__ORB_mdblevec,
               _0RL_MATHEMATICS__ORB_mdblevec_marshal_fn,
               _0RL_MATHEMATICS__ORB_mdblevec_destructor_fn,
               _p);
}
void operator<<=(::CORBA::Any& _a, MATHEMATICS_ORB::dblevec* _sp)
{
  _a.PR_insert(_0RL_tc_MATHEMATICS__ORB_mdblevec,
               _0RL_MATHEMATICS__ORB_mdblevec_marshal_fn,
               _0RL_MATHEMATICS__ORB_mdblevec_destructor_fn,
               _sp);
}

::CORBA::Boolean operator>>=(const ::CORBA::Any& _a, MATHEMATICS_ORB::dblevec*& _sp)
{
  return _a >>= (const MATHEMATICS_ORB::dblevec*&) _sp;
}
::CORBA::Boolean operator>>=(const ::CORBA::Any& _a, const MATHEMATICS_ORB::dblevec*& _sp)
{
  void* _v;
  if (_a.PR_extract(_0RL_tc_MATHEMATICS__ORB_mdblevec,
                    _0RL_MATHEMATICS__ORB_mdblevec_unmarshal_fn,
                    _0RL_MATHEMATICS__ORB_mdblevec_marshal_fn,
                    _0RL_MATHEMATICS__ORB_mdblevec_destructor_fn,
                    _v)) {
    _sp = (const MATHEMATICS_ORB::dblevec*)_v;
    return 1;
  }
  return 0;
}

static void _0RL_MATHEMATICS__ORB_mintvec_marshal_fn(cdrStream& _s, void* _v)
{
  MATHEMATICS_ORB::intvec* _p = (MATHEMATICS_ORB::intvec*)_v;
  *_p >>= _s;
}
static void _0RL_MATHEMATICS__ORB_mintvec_unmarshal_fn(cdrStream& _s, void*& _v)
{
  MATHEMATICS_ORB::intvec* _p = new MATHEMATICS_ORB::intvec;
  *_p <<= _s;
  _v = _p;
}
static void _0RL_MATHEMATICS__ORB_mintvec_destructor_fn(void* _v)
{
  MATHEMATICS_ORB::intvec* _p = (MATHEMATICS_ORB::intvec*)_v;
  delete _p;
}

void operator<<=(::CORBA::Any& _a, const MATHEMATICS_ORB::intvec& _s)
{
  MATHEMATICS_ORB::intvec* _p = new MATHEMATICS_ORB::intvec(_s);
  _a.PR_insert(_0RL_tc_MATHEMATICS__ORB_mintvec,
               _0RL_MATHEMATICS__ORB_mintvec_marshal_fn,
               _0RL_MATHEMATICS__ORB_mintvec_destructor_fn,
               _p);
}
void operator<<=(::CORBA::Any& _a, MATHEMATICS_ORB::intvec* _sp)
{
  _a.PR_insert(_0RL_tc_MATHEMATICS__ORB_mintvec,
               _0RL_MATHEMATICS__ORB_mintvec_marshal_fn,
               _0RL_MATHEMATICS__ORB_mintvec_destructor_fn,
               _sp);
}

::CORBA::Boolean operator>>=(const ::CORBA::Any& _a, MATHEMATICS_ORB::intvec*& _sp)
{
  return _a >>= (const MATHEMATICS_ORB::intvec*&) _sp;
}
::CORBA::Boolean operator>>=(const ::CORBA::Any& _a, const MATHEMATICS_ORB::intvec*& _sp)
{
  void* _v;
  if (_a.PR_extract(_0RL_tc_MATHEMATICS__ORB_mintvec,
                    _0RL_MATHEMATICS__ORB_mintvec_unmarshal_fn,
                    _0RL_MATHEMATICS__ORB_mintvec_marshal_fn,
                    _0RL_MATHEMATICS__ORB_mintvec_destructor_fn,
                    _v)) {
    _sp = (const MATHEMATICS_ORB::intvec*)_v;
    return 1;
  }
  return 0;
}

static void _0RL_MATHEMATICS__ORB_mMATH__Gen_marshal_fn(cdrStream& _s, void* _v)
{
  omniObjRef* _o = (omniObjRef*)_v;
  omniObjRef::_marshal(_o, _s);
}
static void _0RL_MATHEMATICS__ORB_mMATH__Gen_unmarshal_fn(cdrStream& _s, void*& _v)
{
  omniObjRef* _o = omniObjRef::_unMarshal(MATHEMATICS_ORB::MATH_Gen::_PD_repoId, _s);
  _v = _o;
}
static void _0RL_MATHEMATICS__ORB_mMATH__Gen_destructor_fn(void* _v)
{
  omniObjRef* _o = (omniObjRef*)_v;
  if (_o)
    omni::releaseObjRef(_o);
}

void operator<<=(::CORBA::Any& _a, MATHEMATICS_ORB::MATH_Gen_ptr _o)
{
  MATHEMATICS_ORB::MATH_Gen_ptr _no = MATHEMATICS_ORB::MATH_Gen::_duplicate(_o);
  _a.PR_insert(MATHEMATICS_ORB::_tc_MATH_Gen,
               _0RL_MATHEMATICS__ORB_mMATH__Gen_marshal_fn,
               _0RL_MATHEMATICS__ORB_mMATH__Gen_destructor_fn,
               _no->_PR_getobj());
}
void operator<<=(::CORBA::Any& _a, MATHEMATICS_ORB::MATH_Gen_ptr* _op)
{
  _a.PR_insert(MATHEMATICS_ORB::_tc_MATH_Gen,
               _0RL_MATHEMATICS__ORB_mMATH__Gen_marshal_fn,
               _0RL_MATHEMATICS__ORB_mMATH__Gen_destructor_fn,
               (*_op)->_PR_getobj());
  *_op = MATHEMATICS_ORB::MATH_Gen::_nil();
}

::CORBA::Boolean operator>>=(const ::CORBA::Any& _a, MATHEMATICS_ORB::MATH_Gen_ptr& _o)
{
  void* _v;
  if (_a.PR_extract(MATHEMATICS_ORB::_tc_MATH_Gen,
                    _0RL_MATHEMATICS__ORB_mMATH__Gen_unmarshal_fn,
                    _0RL_MATHEMATICS__ORB_mMATH__Gen_marshal_fn,
                    _0RL_MATHEMATICS__ORB_mMATH__Gen_destructor_fn,
                    _v)) {
    omniObjRef* _r = (omniObjRef*)_v;
    if (_r)
      _o = (MATHEMATICS_ORB::MATH_Gen_ptr)_r->_ptrToObjRef(MATHEMATICS_ORB::MATH_Gen::_PD_repoId);
    else
      _o = MATHEMATICS_ORB::MATH_Gen::_nil();
    return 1;
  }
  return 0;
}

