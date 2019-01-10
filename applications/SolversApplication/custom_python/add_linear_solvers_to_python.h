//
//   Project Name:        KratosSolversApplication    $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:         January 2019 $
//   Revision:            $Revision:              0.0 $
//
//

#if !defined(KRATOS_ADD_CUSTOM_LINEAR_SOLVERS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CUSTOM_LINEAR_SOLVERS_TO_PYTHON_H_INCLUDED


// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes


namespace Kratos
{

namespace Python
{

void  AddLinearSolversToPython(pybind11::module& m);

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_CUSTOM_LINEAR_SOLVERS_TO_PYTHON_H_INCLUDED  defined
