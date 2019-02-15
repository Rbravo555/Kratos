//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED )
#define  KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_models/elasticity_models/isochoric_neo_hookean_lnJ_squared_model.hpp"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
  ///@{

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

  /// Short class definition.
  /** Detail class definition.
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IncompressibleNeoHookeanLnJSquaredModel : public IsochoricNeoHookeanLnJSquaredModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IncompressibleNeoHookeanLnJSquaredModel
    KRATOS_CLASS_POINTER_DEFINITION( IncompressibleNeoHookeanLnJSquaredModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IncompressibleNeoHookeanLnJSquaredModel() : IsochoricNeoHookeanLnJSquaredModel() {}

    /// Copy constructor.
    IncompressibleNeoHookeanLnJSquaredModel(IncompressibleNeoHookeanLnJSquaredModel const& rOther) : IsochoricNeoHookeanLnJSquaredModel(rOther) {}

    /// Assignment operator.
    IncompressibleNeoHookeanLnJSquaredModel& operator=(IncompressibleNeoHookeanLnJSquaredModel const& rOther)
    {
      IsochoricNeoHookeanLnJSquaredModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<IncompressibleNeoHookeanLnJSquaredModel>(*this);
    }

    /// Destructor.
    ~IncompressibleNeoHookeanLnJSquaredModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    // Simplyfied methods must be implemented for performance purposes
    /**
     * Calculate Stresses
     */

    /**
     * Calculate Constitutive Components
     */

    /**
     * Check
     */
    int Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      IsochoricNeoHookeanLnJSquaredModel::Check(rProperties,rCurrentProcessInfo);

      return 0;

      KRATOS_CATCH(" ")
    };


    ///@}
    ///@name Access
    ///@{

    /**
     * method to ask the constitutive model the list of variables (dofs) needed from the domain
     * @param rScalarVariables : list of scalar dofs
     * @param rComponentVariables :  list of vector dofs
     */
    void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
                                std::vector<Variable<array_1d<double,3> > >& rComponentVariables) override
    {
      KRATOS_TRY

      HyperElasticModel::GetDomainVariablesList(rScalarVariables, rComponentVariables);

      rScalarVariables.push_back(PRESSURE);

      KRATOS_CATCH(" ")
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IncompressibleHyperElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressibleHyperElasticModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IncompressibleHyperElasticModel Data";
    }


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

    //specialized methods:

    void CalculateVolumetricFactor(HyperElasticDataType& rVariables, double& rFactor) override
    {
      KRATOS_TRY

      rFactor = 1.0;

      KRATOS_CATCH(" ")
    }

    void CalculatePressureFactor(HyperElasticDataType& rVariables, double& rFactor) override
    {
      KRATOS_TRY

      this->CalculateVolumetricFactor(rVariables,rFactor);

      rFactor *= rVariables.GetModelData().GetPressure() * rVariables.Strain.Invariants.J;

      KRATOS_CATCH(" ")
    }

    void CalculateConstitutiveMatrixFactor(HyperElasticDataType& rVariables, double& rFactor) override
    {
      KRATOS_TRY

      rFactor = 1.0;

      KRATOS_CATCH(" ")
    }

    void CalculateConstitutiveMatrixPressureFactor(HyperElasticDataType& rVariables, double& rFactor) override
    {
      KRATOS_TRY

      rFactor = rVariables.GetModelData().GetPressure() * rVariables.Strain.Invariants.J;

      KRATOS_CATCH(" ")
    }


    //************// dW

    double& GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //dU/dJ
    {
      KRATOS_TRY

      const ModelDataType&  rValues = rVariables.GetModelData();

      rDerivative = rValues.GetPressure();

      return rDerivative;

      KRATOS_CATCH(" ")
    };


    double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddU/dJdJ
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };


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
    ///@name Serialization
    ///@{
    friend class Serializer;


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsochoricNeoHookeanLnJSquaredModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsochoricNeoHookeanLnJSquaredModel )
    }


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IncompressibleNeoHookeanLnJSquaredModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_LNJ_SQUARED_MODEL_H_INCLUDED  defined
