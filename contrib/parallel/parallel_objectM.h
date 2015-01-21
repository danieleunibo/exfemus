#ifndef LIBMESH_PARALLEL_OBJECT_H
#define LIBMESH_PARALLEL_OBJECT_H




// Local includes
#include "parallelM.h"

// mydefine --------------------------------------------------------
#define processor_id_type int
// Macro to identify and debug functions which should only be called in
// parallel on every processor at once
#undef parallel_object_onlyM
#ifndef NDEBUG
#define parallel_object_onlyM() parallel_onlyM(this->comm())
#else
#define parallel_object_onlyM()  ((void) 0)
#endif
// -----------------------------------------------------------------

// =================================================================
///  This class forms the base class for all other classes
///  that are expected to be implemented in paralel. Each
///  \p ParalelObject *requires* a \p ParallelM::Communicator object
///  for construction.
class ParallelObjectM {
// =================================================================
protected:
    // data ------------------------------------------
    const ParallelM::Communicator &_communicator;

public:
    // Constructor- Destructor ----------------------
    /// Constructor. Requires a reference to the communicator
    ParallelObjectM (
        const ParallelM::Communicator &comm_in
    ) :  _communicator(comm_in) {}

    //// Copy Constructor.
    ParallelObjectM (
        const ParallelObjectM &other
    ) :  _communicator(other._communicator) {}

    /// "Assignment" operator.  Simply asserts our references
    ParallelObjectM & operator= (
        const ParallelObjectM &other
    ) {
        assert(&_communicator == &other._communicator);
        return *this;
    }

    /// Destructor.  Virtual because we are a base class.
    virtual ~ParallelObjectM () {}

    // Return functions --------------------------------
    // @returns a reference to the \p ParallelM::Communicator object used by this mesh.
    const ParallelM::Communicator & comm() const    {
        return _communicator;
    }
    ///  @returns the number of processors in the group.
    processor_id_type n_processors () const {
        return static_cast<processor_id_type>(_communicator.size());
    }
    /// @returns the rank of this processor in the group.
    processor_id_type processor_id () const {
        return static_cast<processor_id_type>(_communicator.rank());
    }

};

#endif // LIBMESH_PARALLEL_OBJECT_H
