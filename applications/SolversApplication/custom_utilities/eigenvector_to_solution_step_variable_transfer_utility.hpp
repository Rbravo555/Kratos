//
//   Project Name:        KratosSolversApplication $
//   Developed by:        $Developer:  JMCarbonell $
//   Maintained by:       $Maintainer:        JMC  $
//   Date:                $Date:      January 2019 $
//
//

#if !defined(KRATOS_EIGENVECTOR_TO_SOLUTION_STEP_VARIABLE_TRANSFER_UTILITY_H_INCLUDED)
#define  KRATOS_EIGENVECTOR_TO_SOLUTION_STEP_VARIABLE_TRANSFER_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes
#include "solvers_application_variables.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Transfer eigenvectors to solution step variables for GiD output or solution initialization.
/**
 * Example Python Code:
 * # Eigenvectors are first computed and stored in the nodal variable EIGENVECTOR_MATRIX.
 * for step in range(NumEigenvalues):
 *   main_model_part.ProcessInfo[TIME] = float(step+1)
 *   EigenvectorToSolutionStepVariableTransferUtility().Transfer(main_model_part,step,0)
 *   gid_output.PrintOutput()
 */
class EigenvectorToSolutionStepVariableTransferUtility {
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( EigenvectorToSolutionStepVariableTransferUtility );

    ///@}
    ///@name Life Cycle
    ///@{

    EigenvectorToSolutionStepVariableTransferUtility() {}

    virtual ~EigenvectorToSolutionStepVariableTransferUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Transfer(ModelPart& rModelPart, int iEigenMode, int step=0)
    {
        for (auto itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++)
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = itNode->GetDofs();
            Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);
            std::size_t j=0;
            for (auto itDof = std::begin(rNodeDofs); itDof != std::end(rNodeDofs); itDof++)
                itDof->GetSolutionStepValue(step) = rNodeEigenvectors(iEigenMode,j++);
        }
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // class EigenvectorToSolutionStepVariableTransferUtility

///@}

///@name Type Definitions
///@{

///@}

}
 // namespace Kratos
#endif  // KRATOS_EIGENVECTOR_TO_SOLUTION_STEP_VARIABLE_TRANSFER_UTILITY_H_INCLUDED defined
