//
//   Project Name:        KratosSolversApplication $
//   Developed by:        $Developer:  JMCarbonell $
//   Maintained by:       $Maintainer:        JMC  $
//   Date:                $Date:      January 2019 $
//
//

#if !defined(KRATOS_NEWMARK_STEP_ROTATION_METHOD_H_INCLUDED)
#define  KRATOS_NEWMARK_STEP_ROTATION_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/time_integration_methods/newmark_step_method.hpp"

namespace Kratos
{
  ///@addtogroup SolversApplication
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
   * This class performs predict and update of dofs variables, their time derivatives and time integrals
   */
  template<class TVariableType, class TValueType>
  class KRATOS_API(SOLVERS_APPLICATION) NewmarkStepRotationMethod : public NewmarkStepMethod<TVariableType,TValueType>
  {
  public:

    ///@name Type Definitions
    ///@{

    /// BaseType
    typedef TimeIntegrationMethod<TVariableType,TValueType>  BaseType;

    /// BasePointerType
    typedef typename BaseType::Pointer                BasePointerType;

    /// NodeType
    typedef typename BaseType::NodeType                      NodeType;

    /// KratosVariable or KratosVariableComponent
    typedef typename BaseType::VariablePointer        VariablePointer;

    /// DerivedType
    typedef NewmarkStepMethod<TVariableType,TValueType>   DerivedType;


    KRATOS_CLASS_POINTER_DEFINITION( NewmarkStepRotationMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    NewmarkStepRotationMethod() : DerivedType() {}

    /// Constructor.
    NewmarkStepRotationMethod(const TVariableType& rVariable) : DerivedType(rVariable) {}

    /// Constructor.
    NewmarkStepRotationMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative) {}

    /// Constructor.
    NewmarkStepRotationMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable) {}


    /// Copy Constructor.
    NewmarkStepRotationMethod(NewmarkStepRotationMethod& rOther)
      :DerivedType(rOther)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new NewmarkStepRotationMethod(*this) );
    }

    /// Destructor.
    ~NewmarkStepRotationMethod() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    // update
    void Update(NodeType& rNode) override;

    ///@}
    ///@name Access
    ///@{

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
        buffer << "NewmarkStepRotationMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NewmarkStepRotationMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NewmarkStepRotationMethod Data";
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
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
    };

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )
    };

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NewmarkStepRotationMethod

  ///@}

  ///@name Type Definitions
  ///@{

  template<>
  void NewmarkStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode);

  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, NewmarkStepRotationMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const NewmarkStepRotationMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEWMARK_STEP_ROTATION_METHOD_H_INCLUDED defined
