//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_VON_MISES_NEO_HOOKEAN_LINEAR_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_VON_MISES_NEO_HOOKEAN_LINEAR_PLASTICITY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/linear_associative_plasticity_model.hpp"
#include "custom_models/plasticity_models/yield_surfaces/mises_huber_yield_surface.hpp"
#include "custom_models/plasticity_models/hardening_rules/simo_linear_hardening_rule.hpp"
#include "custom_models/elasticity_models/isochoric_neo_hookean_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) VonMisesNeoHookeanLinearPlasticityModel : public LinearAssociativePlasticityModel<IsochoricNeoHookeanModel, MisesHuberYieldSurface<SimoLinearHardeningRule> >
  {
  public:

    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef IsochoricNeoHookeanModel                       ElasticityModelType;
    typedef ElasticityModelType::Pointer                ElasticityModelPointer;

    //yield surface
    typedef SimoLinearHardeningRule                          HardeningRuleType;
    typedef MisesHuberYieldSurface<HardeningRuleType>         YieldSurfaceType;
    typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

    //base type
    typedef LinearAssociativePlasticityModel<ElasticityModelType,YieldSurfaceType>  BaseType;

    //common types
    typedef BaseType::Pointer                         BaseTypePointer;
    typedef BaseType::SizeType                               SizeType;
    typedef BaseType::VoigtIndexType                   VoigtIndexType;
    typedef BaseType::MatrixType                           MatrixType;
    typedef BaseType::ModelDataType                     ModelDataType;
    typedef BaseType::MaterialDataType               MaterialDataType;
    typedef BaseType::PlasticDataType                 PlasticDataType;
    typedef BaseType::InternalVariablesType     InternalVariablesType;


    /// Pointer definition of VonMisesNeoHookeanLinearPlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( VonMisesNeoHookeanLinearPlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VonMisesNeoHookeanLinearPlasticityModel() : BaseType() {}

    /// Copy constructor.
    VonMisesNeoHookeanLinearPlasticityModel(VonMisesNeoHookeanLinearPlasticityModel const& rOther) : BaseType(rOther) {}

    /// Assignment operator.
    VonMisesNeoHookeanLinearPlasticityModel& operator=(VonMisesNeoHookeanLinearPlasticityModel const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<VonMisesNeoHookeanLinearPlasticityModel>(*this);
    }

    /// Destructor.
    ~VonMisesNeoHookeanLinearPlasticityModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    /**
     * Has Values
     */
    bool Has(const Variable<double>& rThisVariable) override
    {
      if(rThisVariable == PLASTIC_STRAIN || rThisVariable == DELTA_PLASTIC_STRAIN )
	return true;

      return false;
    }


    /**
     * Get Values
     */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
    {

      rValue=0;

      if (rThisVariable==PLASTIC_STRAIN)
	{
	  rValue = this->mInternal.Variables[0];
	}


      if (rThisVariable==DELTA_PLASTIC_STRAIN)
	{
	  rValue = this->mInternal.Variables[0]-mPreviousInternal.Variables[0];
	}


      return rValue;
    }

    /**
     * Set Values
     */
    void SetValue(const Variable<double>& rThisVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override
    {
      if (rThisVariable==PLASTIC_STRAIN)
	{
            this->mInternal.Variables[0] = rValue;
	}
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
      buffer << "VonMisesNeoHookeanLinearPlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesNeoHookeanLinearPlasticityModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesNeoHookeanLinearPlasticityModel Data";
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
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class VonMisesNeoHookeanLinearPlasticityModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_VON_MISES_NEO_HOOKEAN_LINEAR_PLASTICITY_MODEL_H_INCLUDED  defined
