//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              LMonforte $
//   Last modified by:    $Co-Author:         JMCarbonell $
//   Date:                $Date:             January 2019 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_J_ELEMENT_H_INCLUDED )
#define KRATOS_UPDATED_LAGRANGIAN_U_J_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/large_displacement_element.hpp"

namespace Kratos
{
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

/// Updated Lagrangian U-J Element for 3D and 2D geometries. Linear Triangles and Tetrahedra

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for Linear Triangles and Tetrahedra.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) UpdatedLagrangianUJElement
    : public LargeDisplacementElement
{
public:

  ///@name Type Definitions
  ///@{
  ///Reference type definition for constitutive laws
  typedef ConstitutiveLaw ConstitutiveLawType;
  ///Pointer type for constitutive laws
  typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
  ///StressMeasure from constitutive laws
  typedef ConstitutiveLawType::StressMeasure StressMeasureType;
  ///Type definition for integration methods
  typedef GeometryData::IntegrationMethod IntegrationMethod;
  ///Type for size
  typedef GeometryData::SizeType SizeType;
  ///Type for element variables
  typedef LargeDisplacementElement::ElementDataType ElementDataType;

  /// Counted pointer of LargeDisplacementUPElement
  KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUJElement );

  ///@}
  ///@name Life Cycle
  ///@{


  /// Default constructors
  UpdatedLagrangianUJElement(IndexType NewId, GeometryType::Pointer pGeometry);

  UpdatedLagrangianUJElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

  ///Copy constructor
  UpdatedLagrangianUJElement(UpdatedLagrangianUJElement const& rOther);

  /// Destructor.
  ~UpdatedLagrangianUJElement() override;

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  UpdatedLagrangianUJElement& operator=(UpdatedLagrangianUJElement const& rOther);

  ///@}
  ///@name Operations
  ///@{

  /**
   * creates a new total lagrangian updated element pointer
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

  /**
   * clones the selected element variables, creating a new one
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

  //************* GETTING METHODS

  //SET

  /**
   * Set a double  Value on the Element Constitutive Law
   */
  void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  //GET:

  /**
   * Get on rVariable a double Value from the Element Constitutive Law
   */
  void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo) override;


  //on integration points:

  /**
   * Calculate a double Variable on the Element Constitutive Law
   */
  void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

  void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

  void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;



  //************* STARTING - ENDING  METHODS

  /**
   * Called to initialize the element.
   * Must be called before any calculation is done
   */
  void Initialize() override;

  /**
   * Sets on rElementalDofList the degrees of freedom of the considered element geometry
   */
  void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Sets on rResult the ID's of the element degrees of freedom
   */
  void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Sets on rValues the nodal displacements
   */
  void GetValuesVector(Vector& rValues, int Step = 0) override;

  /**
   * Sets on rValues the nodal velocities
   */
  void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

  /**
   * Sets on rValues the nodal accelerations
   */
  void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

  //************* COMPUTING  METHODS

  /**
   * this is called during the assembling process in order
   * to calculate the elemental mass matrix
   * @param rMassMatrix: the elemental mass matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the elemental damping matrix
   * @param rDampingMatrix: the elemental damping matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;

  //************************************************************************************
  //************************************************************************************
  /**
   * This function provides the place to perform checks on the completeness of the input.
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   */
  int Check(const ProcessInfo& rCurrentProcessInfo) override;

  ///@}
  ///@name Access
  ///@{
  ///@}
  ///@name Inquiry
  ///@{
  ///@}
  ///@name Input and output
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

  /**
   * Container for historical total elastic deformation measure F0 = dx/dX
   */
  std::vector< Matrix > mDeformationGradientF0;
  std::vector< Matrix > mDeformationGradientJ0;

  /**
   * Container for the total deformation gradient determinants
   */
  Vector mDeterminantF0;
  Vector mDeterminantJ0;

  ///@}
  ///@name Protected Operators
  ///@{

  /// Empty constructor needed for serialization
  UpdatedLagrangianUJElement();

  ///@}
  ///@name Protected Operations
  ///@{

  /**
   * Get dof size of a node
   */
  SizeType GetNodeDofsSize() override;

  /**
   * Initialize Element General Variables
   */
  void InitializeElementData(ElementDataType & rVariables,
                             const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Finalize Element Internal Variables
   */
  void FinalizeStepVariables(ElementDataType & rVariables,
                             const double& rPointNumber) override;

  /**
   * Calculate Element Kinematics
   */
  void CalculateKinematics(ElementDataType& rVariables,
                           const double& rPointNumber) override;

  /**
   * Calculation and addition of the matrices of the LHS
   */
  void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                          ElementDataType& rVariables,
                          double& rIntegrationWeight) override;

  /**
   * Calculation and addition of the vectors of the RHS
   */
  void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                          ElementDataType& rVariables,
                          Vector& rVolumeForce,
                          double& rIntegrationWeight) override;

  /**
   * Calculation of the Internal Forces due to sigma. Fi = B * sigma
   */
  void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
                                     ElementDataType & rVariables,
                                     double& rIntegrationWeight) override;

  /**
   * Calculation of the Internal Forces due to Pressure-Balance
   */
  virtual void CalculateAndAddJacobianForces(VectorType& rRightHandSideVector,
                                             ElementDataType & rVariables,
                                             double& rIntegrationWeight);

  /**
   * Calculation of the Internal Forces due to Pressure-Balance
   */
  virtual void CalculateAndAddStabilizedJacobian(VectorType& rRightHandSideVector,
                                                 ElementDataType & rVariables,
                                                 double& rIntegrationWeight);

  /**
   * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
   */
  void CalculateAndAddKuum(MatrixType& rK,
                           ElementDataType & rVariables,
                           double& rIntegrationWeight) override;

  /**
   * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
   */
  void CalculateAndAddKuug(MatrixType& rK,
                           ElementDataType & rVariables,
                           double& rIntegrationWeight) override;

  /**
   * Calculation of the Kup matrix
   */
  virtual void CalculateAndAddKuJ(MatrixType& rK,
                                  ElementDataType & rVariables,
                                  double& rIntegrationWeight);

  /**
   * Calculation of the Kpu matrix
   */
  virtual void CalculateAndAddKJu(MatrixType& rK,
                                  ElementDataType & rVariables,
                                  double& rIntegrationWeight);

  /**
   * Calculation of the Kpp matrix
   */
  virtual void CalculateAndAddKJJ(MatrixType& rK,
                                  ElementDataType & rVariables,
                                  double& rIntegrationWeight);

  /**
   * Calculation of the Kpp Stabilization Term matrix
   */
  virtual void CalculateAndAddKJJStab(MatrixType& rK, ElementDataType & rVariables,
                                      double& rIntegrationWeight);

  /**
   * Get the Historical Deformation Gradient to calculate after finalize the step
   */
  void GetHistoricalVariables(ElementDataType& rVariables,
                              const double& rPointNumber) override;

  /**
   * Get Stabilization Parameter
   */
  virtual void CalculateStabilizationParameter(double& rStabilizationParameter, const ProcessInfo& rCurrentProcessInfo);

  /**
   * Calculation of the Volume Change of the Element
   */
  double& CalculateVolumeChange(double& rVolumeChange,
                                ElementDataType& rVariables) override;


  /**
   * Set Variables of the Element to the Parameters of the Constitutive Law
   */
  void SetElementData(ElementDataType& rVariables,
                      ConstitutiveLaw::Parameters& rValues,
                      const int & rPointNumber) override;
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
  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  // A private default constructor necessary for serialization

  void save(Serializer& rSerializer) const override;

  void load(Serializer& rSerializer) override;

  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}


}; // Class UpdatedLagrangianUJElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_ELEMENT_H_INCLUDED
