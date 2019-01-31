//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              LMonforte $
//   Last modified by:    $Co-Author:         JMCarbonell $
//   Date:                $Date:             January 2019 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_J_element.hpp"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUJElement::UpdatedLagrangianUJElement()
    : LargeDisplacementElement()
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUJElement::UpdatedLagrangianUJElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUJElement::UpdatedLagrangianUJElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianUJElement::UpdatedLagrangianUJElement( UpdatedLagrangianUJElement const& rOther)
    :LargeDisplacementElement(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeformationGradientJ0(rOther.mDeformationGradientJ0)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mDeterminantJ0(rOther.mDeterminantJ0)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianUJElement&  UpdatedLagrangianUJElement::operator=(UpdatedLagrangianUJElement const& rOther)
{
  LargeDisplacementElement::operator=(rOther);

  mDeformationGradientF0.clear();
  mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

  mDeformationGradientJ0.clear();
  mDeformationGradientJ0.resize(rOther.mDeformationGradientJ0.size());

  for(SizeType i=0; i<mConstitutiveLawVector.size(); ++i)
  {
    mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
    mDeformationGradientJ0[i] = rOther.mDeformationGradientJ0[i];
  }

  mDeterminantF0 = rOther.mDeterminantF0;
  mDeterminantJ0 = rOther.mDeterminantJ0;

  return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUJElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared< UpdatedLagrangianUJElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUJElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

  Element::Pointer pClonedElement = Kratos::make_shared<UpdatedLagrangianUJElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
  //-----------//

  pClonedElement.mThisIntegrationMethod = mThisIntegrationMethod;


  if ( pClonedElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
  {
    pClonedElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

    if( pClonedElement.mConstitutiveLawVector.size() != pClonedElement.GetGeometry().IntegrationPointsNumber() )
      KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", pClonedElement.mConstitutiveLawVector.size() )
          }

  for(SizeType i=0; i<mConstitutiveLawVector.size(); i++)
  {
    pClonedElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
  }

  //-----------//

  if ( pClonedElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
    pClonedElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

  if ( pClonedElement.mDeformationGradientJ0.size() != mDeformationGradientJ0.size() )
    pClonedElement.mDeformationGradientJ0.resize(mDeformationGradientJ0.size());


  for(SizeType i=0; i<mDeformationGradientJ0.size(); i++)
  {
    pClonedElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
    pClonedElement.mDeformationGradientJ[i] = mDeformationGradientJ0[i];
  }

  pClonedElement.mDeterminantF0 = mDeterminantF0;
  pClonedElement.mDeterminantJ0 = mDeterminantJ0;

  pClonedElement.SetData(this->GetData());
  pClonedElement.SetFlags(this->GetFlags());

  return pClonedElement;
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUJElement::~UpdatedLagrangianUJElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
  rElementalDofList.resize( 0 );

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  for ( SizeType i = 0; i < GetGeometry().size(); i++ )
  {
    rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
    rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

    if( dimension == 3 )
      rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

    rElementalDofList.push_back( GetGeometry()[i].pGetDof( JACOBIAN ));
  }
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType dofs_size = this->GetDofsSize();

  if ( rResult.size() != dofs_size )
    rResult.resize( dofs_size, false );

  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    int index = i * dimension + i;
    rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
    rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

    if( dimension == 3)
    {
      rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
      rResult[index + 3] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
    }
    else
    {
      rResult[index + 2] = GetGeometry()[i].GetDof( JACOBIAN ).EquationId();
    }

  }

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void LargeDisplacementUJElement::GetValuesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = number_of_nodes * dimension + number_of_nodes;

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        SizeType index = i * dimension + i;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
        }
        else
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( JACOBIAN, Step );
        }

    }
}


//************************************VELOCITY****************************************
//************************************************************************************

void LargeDisplacementUJElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = number_of_nodes * dimension + number_of_nodes;

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        SizeType index = i * dimension + i;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
        if ( dimension == 3 )
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index + 3] = 0;
        }
        else
        {
            rValues[index + 2] = 0;
        }
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void LargeDisplacementUJElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = number_of_nodes * dimension + number_of_nodes;

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        SizeType index = i * dimension + i;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index + 3] = 0;
        }
        else
        {
            rValues[index + 2] = 0;
        }
    }

}

//**************************************************************************
//**************************************************************************

UpdatedLagrangianUJElement::SizeType UpdatedLagrangianUJElement::GetNodeDofsSize()
{
  return (GetGeometry().WorkingSpaceDimension() + 1); //usual size for U-J elements
}

//**************************************************************************
//**************************************************************************

void UpdatedLagrangianUJElement::CalculateStabilizationFactor(double& rStabilizationFactor)
{
  KRATOS_TRY

  //LMonforte
  double YoungModulus = 1.0;
  if( GetProperties().Has(YOUNG_MODULUS) )
    YoungModulus = GetProperties()[YOUNG_MODULUS];

  if(YoungModulus<1e-5)
  {
    double BulkModulus = 0.0;
    if(GetProperties().Has(BULK_MODULUS)){
      BulkModulus= GetProperties()[BULK_MODULUS];
    }
    else if( GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO) ){
      BulkModulus = GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
    }

    double ShearModulus = 0.0;
    double BulkModulus = 1.0;
    if(GetProperties().Has(SHEAR_MODULUS)){
      ShearModulus= GetProperties()[SHEAR_MODULUS];
    }
    else if(GetProperties().Has(BULK_MODULUS) && GetProperties().Has(POISSON_RATIO)){
      ShearModulus= GetProperties()[BULK_MODULUS]*3*(1-2*GetProperties()[POISSON_RATIO])/(2*(1+GetProperties()[POISSON_RATIO]));
    }
    else if( GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO) ){
      ShearModulus = GetProperties()[YOUNG_MODULUS]/(2*(1+GetProperties()[POISSON_RATIO]));
    }

    if(fabs(ShearModulus) > 1e-6)
      rStabilizationFactor = 1.0 / ShearModulus;

    GetValueOnIntegrationPoints(BULK_MODULUS, Values, rCurrentProcessInfo);
    if(fabs(BulkModulus) > 1e-6)
      rStabilizationFactor *= BulkModulus;
  }

  KRATOS_CATCH("");
}


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************
void UpdatedLagrangianUJElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                              std::vector<double>& rValues,
                                                              const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const SizeType& integration_points_number = mConstitutiveLawVector.size();
    for ( SizeType PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
    {
      mDeterminantF0[PointNumber] = rValues[PointNumber];
    }

  }
  else{

    LargeDisplacementElement::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }


}


//**********************************GET DOUBLE VALUE**********************************
//************************************************************************************
void UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                              std::vector<double>& rValues,
                                                              const ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  const SizeType& integration_points_number = mConstitutiveLawVector.size();
  if (rVariable == DETERMINANT_F){

    if ( rValues.size() != integration_points_number )
      rValues.resize( integration_points_number );

    for ( SizeType PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
    {
      rValues[PointNumber] = mDeterminantF0[PointNumber];
    }
  }
  else{
    LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
  }

  KRATOS_CATCH("")
}

// ******************************************************************************************
// Calculate On Integration Points. (VECTOR)
void UpdatedLagrangianUJElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{

  KRATOS_TRY

      const SizeType& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

  if ( rOutput.size() != integration_points_number )
    rOutput.resize( integration_points_number );
  if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
  {
    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

    //reading integration points
    for ( SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
      //compute element kinematics B, F, DN_DX ...
      this->CalculateKinematics(Variables,PointNumber);

      //to take in account previous step writing
      if( this->Is(SolidElement::FINALIZED_STEP) ){
        this->GetHistoricalVariables(Variables,PointNumber);
      }

      //set general variables to constitutivelaw parameters
      this->SetElementData(Variables,Values,PointNumber);

      // OBS, now changing Variables I change Values because they are pointers ( I hope);
      double ElementalDetFT = Variables.detH;
      Matrix ElementalFT = Variables.H;

      // AND NOW IN THE OTHER WAY
      Matrix m; double d;
      this->ComputeConstitutiveVariables( Variables, m, d);

      Variables.H = m;
      Variables.detH = d;
      Values.SetDeformationGradientF( Variables.H);
      Values.SetDeterminantF( Variables.detH );

      //call the constitutive law to update material variables
      if( rVariable == CAUCHY_STRESS_VECTOR)
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
      else
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);

      Variables.H = ElementalFT;
      Variables.detH = ElementalDetFT;

      if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
        rOutput[PointNumber].resize( Variables.StressVector.size(), false );

      rOutput[PointNumber] = Variables.StressVector;

    }

  }
  else if ( rVariable == DARCY_FLOW)
  {
    Properties thisProperties = GetProperties();
    // CONSTITUTIVE PARAMETERS
    double Permeability = 0;
    if( GetProperties().Has(PERMEABILITY) ){
      Permeability = GetProperties()[PERMEABILITY];
    }
    double WaterDensity = 0;
    if( GetProperties().Has(DENSITY_WATER) ){
      WaterDensity = GetProperties()[DENSITY_WATER];
    }

    // GEOMETRY PARAMETERS
    const SizeType& integration_points_number = mConstitutiveLawVector.size();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();
    const SizeType& number_of_nodes = GetGeometry().size();

    // Get DN_DX
    ElementDataType Variables;
    this->InitializeElementData( Variables, rCurrentProcessInfo);

    Matrix K = ZeroMatrix( dimension, dimension);
    for (SizeType i = 0; i < dimension; i++)
      K(i,i) = Permeability;  // this is only one of the two cases.

    for (SizeType PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
    {
      this->CalculateKinematics(Variables, PointNumber);

      Vector GradP = ZeroVector( dimension );

      for (SizeType i = 0; i < number_of_nodes; i++) {
        if ( GetGeometry()[i].HasDofFor( WATER_PRESSURE ) == false) {
          return;
        }
        const double & rWaterPressure = GetGeometry()[i].FastGetSolutionStepValue( WATER_PRESSURE );
        for (SizeType iDim = 0; iDim < dimension; iDim++) {
          GradP(iDim) += Variables.DN_DX(i, iDim) * rWaterPressure;
        }
      }

      // BTerm
      GradP(dimension-1) -= 10.0 * WaterDensity;

      // finally
      GradP  = prod( K, GradP);
      // Manual resize
      Vector ResizedVector = ZeroVector(3);
      for (SizeType i = 0; i < dimension; i++) {
        ResizedVector(i) = GradP(i);
      }
      rOutput[PointNumber] = ResizedVector;

    }
  }
  else {
    LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
  }

  KRATOS_CATCH("")
      }

// ******************************************************************************************
// Calculate On Integration Points. (MATRIX)
void UpdatedLagrangianUJElement::CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

      if (rVariable == CAUCHY_STRESS_TENSOR)
      {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        Variables.StressVector = ZeroVector(6); // I WANT TO GET THE THIRD COMPONENT
        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        //reading integration points
        for ( SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
          //compute element kinematics B, F, DN_DX ...
          this->CalculateKinematics(Variables,PointNumber);

          //to take in account previous step writing
          if( this->Is(SolidElement::FINALIZED_STEP) ){
            this->GetHistoricalVariables(Variables,PointNumber);
          }

          //set general variables to constitutivelaw parameters
          this->SetElementData(Variables,Values,PointNumber);

          // OBS, now changing Variables I change Values because they are pointers ( I hope);
          double ElementalDetFT = Variables.detH;
          Matrix ElementalFT = Variables.H;

          // AND NOW IN THE OTHER WAY
          Matrix m; double d;
          this->ComputeConstitutiveVariables( Variables, m, d);

          Variables.H = m;
          Variables.detH = d;
          Values.SetDeformationGradientF( Variables.H);
          Values.SetDeterminantF( Variables.detH );

          mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

          Variables.H = ElementalFT;
          Variables.detH = ElementalDetFT;

          if ( ( rOutput[PointNumber].size1() != 3 ) |
               ( rOutput[PointNumber].size2() != 3 ) )
            rOutput[PointNumber].resize( 3, 3, false );

          rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor( Variables.StressVector );


        }

      }
      else if ( rVariable == TOTAL_CAUCHY_STRESS) {

        CalculateOnIntegrationPoints( CAUCHY_STRESS_TENSOR, rOutput, rCurrentProcessInfo);

        if ( GetGeometry()[0].HasDofFor( WATER_PRESSURE) )
        {
          const SizeType number_of_nodes = GetGeometry().size();
          //create and initialize element variables:
          ElementDataType Variables;
          this->InitializeElementData(Variables,rCurrentProcessInfo);

          //reading integration points
          for ( SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
          {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            double WaterPressure = 0;
            for (SizeType i = 0; i < number_of_nodes; i++)
            {
              WaterPressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE) * Variables.N[i];
            }

            for (SizeType i = 0; i < 3; i++)
              rOutput[PointNumber](i,i) += WaterPressure;

          }
        }

      }
      else {
        LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }

  KRATOS_CATCH("")
      }

// ******************************************************************************************
// Calculate On Integration Points. (DOUBLE)
void UpdatedLagrangianUJElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

      if ( rVariable == DETERMINANT_F) {
        const SizeType& integration_points_number = mConstitutiveLawVector.size();

        if (rOutput.size() != integration_points_number)
          rOutput.resize( integration_points_number) ;

        for ( SizeType PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
          rOutput[PointNumber] = mDeterminantF0[PointNumber];
        }

      }
      else if ( rVariable == POROSITY)
      {

        const SizeType& integration_points_number = mConstitutiveLawVector.size();
        const double InitialPorosity  = GetProperties()[INITIAL_POROSITY];

        if ( rOutput.size() != mConstitutiveLawVector.size() )
          rOutput.resize( mConstitutiveLawVector.size() );

        std::vector<double>  DetF0;
        GetValueOnIntegrationPoints( DETERMINANT_F, DetF0, rCurrentProcessInfo);

        if (rOutput.size() != integration_points_number)
          rOutput.resize( integration_points_number) ;

        for ( SizeType PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
          rOutput[PointNumber] = 1.0 - (1.0 - InitialPorosity) / DetF0[PointNumber] ;
        }
      }
      else if ( rVariable == VOID_RATIO) {

        GetValueOnIntegrationPoints( POROSITY, rOutput, rCurrentProcessInfo);

        const SizeType& integration_points_number = mConstitutiveLawVector.size();

        for ( SizeType PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
          rOutput[PointNumber] = rOutput[PointNumber] / (1.0 - rOutput[PointNumber]) ;
        }

      }
      else if ( rVariable == PERMEABILITY)
      {
        const SizeType& integration_points_number = mConstitutiveLawVector.size();

        if ( rOutput.size() != mConstitutiveLawVector.size() )
          rOutput.resize( mConstitutiveLawVector.size() );

        double Permeability    = GetProperties()[PERMEABILITY];
        bool Kozeny = GetProperties()[KOZENY_CARMAN];
        if ( Kozeny == false)
        {
          for ( SizeType PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
          {
            rOutput[PointNumber] = Permeability ;
          }
          return;
        }

        std::vector<double>  Porosity;
        GetValueOnIntegrationPoints( POROSITY, Porosity, rCurrentProcessInfo);
        double PorosityInitial = GetProperties()[INITIAL_POROSITY];
        double initialVoidRatio = PorosityInitial / (1.0 - PorosityInitial);

        double Constant = Permeability * ( 1.0 + initialVoidRatio) / pow( initialVoidRatio, 3.0);

        for ( SizeType PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
          double voidRatio = Porosity[PointNumber] / ( 1.0 - Porosity[PointNumber]);
          rOutput[PointNumber] = Constant * pow( voidRatio, 3.0) / (1.0 + voidRatio);
          if ( rOutput[PointNumber] < Permeability / 1000.0) {
            rOutput[PointNumber] = Permeability /1000.0;
          }
        }


      }
      else {
        LargeDisplacementElement::CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo);
      }

  KRATOS_CATCH("")
      }



void UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( const Variable<Vector> & rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
{

  if ( rVariable == DARCY_FLOW)
  {
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo);
  } else {
    LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
  }
}

//**********************************GET TENSOR VALUE**********************************
//************************************************************************************
void UpdatedLagrangianUJElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo)
{
  if ( rVariable == CAUCHY_STRESS_TENSOR) {
    CalculateOnIntegrationPoints(rVariable, rValue, rCurrentProcessInfo);
  }
  else if ( rVariable == TOTAL_CAUCHY_STRESS) {
    CalculateOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);
  }
  else {
    LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValue, rCurrentProcessInfo);
  }

}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************
void UpdatedLagrangianUJElement::Initialize()
{
  KRATOS_TRY

  LargeDisplacementElement::Initialize();

  const SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  //Resize historic deformation gradient
  if(mDeformationGradientF0.size() != integration_points_number)
    mDeformationGradientF0.resize(integration_points_number, false);

  if(mDeterminantF0.size() != integration_points_number)
    mDeterminantF0.resize(integration_points_number, false);

  if(mDeformationGradientJ0.size() != integration_points_number)
    mDeformationGradientJ0.resize(integration_points_number, false);

  if(mDeterminantJ0.size() != integration_points_number)
    mDeterminantJ0.resize(integration_points_number, false);

  for(SizeType PointNumber=0; PointNumber<integration_points_number; ++PointNumber)
  {
    mDeterminantF0[PointNumber] = 1;
    noalias(mDeformationGradientF0[PointNumber]) = IdentityMatrix(dimension);
    mDeterminantJ0[PointNumber] = 1;
    noalias(mDeformationGradientJ0[PointNumber]) = IdentityMatrix(dimension);
  }

  //ATTENTION initialize nodal variables (parallelism)
  const SizeType number_of_nodes = GetGeometry().size();
  for(SizeType i=0; i<number_of_nodes; ++i) {
    GetGeometry()[i].SetLock();
    GetGeometry()[i].FastGetSolutionStepValue(JACOBIAN) = 1.0
    GetGeometry()[i].FastGetSolutionStepValue(JACOBIAN,1) = 1.0;
    GetGeometry()[i].UnSetLock();
  }

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::InitializeElementData (ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  LargeDisplacementElement::InitializeElementData(rVariables,rCurrentProcessInfo);

  //Calculate Delta Position
  ElementUtilities::CalculateDeltaPosition(rVariables.DeltaPosition,this->GetGeometry());

  //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
  rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::FinalizeStepVariables( ElementDataType & rVariables, const double& rPointNumber )
{
  KRATOS_TRY

  //update internal (historical) variables
  mDeterminantF0[rPointNumber] = rVariables.detF * rVariables.detF0;
  noalias(mDeformationGradientF0[rPointNumber]) = prod(rVariables.F, rVariables.F0);

  mDeterminantJ0[rPointNumber] = rVariables.detH;
  noalias(mDeformationGradientJ0[rPointNumber]) = rVariables.H;

  KRATOS_CATCH("")
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************
void UpdatedLagrangianUJElement::CalculateKinematics(ElementDataType& rVariables,
                                                     const double& rPointNumber)
{
  KRATOS_TRY

  //Get the parent coodinates derivative [dN/d£]
  const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

  //Get the shape functions for the order of the integration method [N]
  const Matrix& Ncontainer = rVariables.GetShapeFunctions();

  //Parent to reference configuration
  rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

  //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
  Matrix InvJ;
  MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

  //Compute cartesian derivatives [dN/dx_n]
  noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

  //Deformation Gradient F [dx_n+1/dx_n] to be updated
  noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

  //Determinant of the deformation gradient F
  rVariables.detF = MathUtils<double>::Det(rVariables.F);

  //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
  Matrix Invj;
  MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

  //Compute cartesian derivatives [dN/dx_n+1]
  noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

  //Determinant of the Deformation Gradient F0
  rVariables.detF0 = mDeterminantF0[rPointNumber];
  rVariables.F0    = mDeformationGradientF0[rPointNumber];

  //Set Shape Functions Values for this integration point
  noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

  //Compute the deformation matrix B
  const GeometryType& rGeometry = GetGeometry();
  ElementUtilities::CalculateLinearDeformationMatrix(rVariables.B,rGeometry,rVariables.DN_DX);

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
{
  rVariables.detF0 *= rVariables.detF;
  double DeterminantF = rVariables.detF;
  rVariables.detF = 1.0;

  //contributions of the stiffness matrix calculated on the reference configuration
  MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

  // operation performed: add Km to the rLefsHandSideMatrix

  //respect to the current configuration n+1
  this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  // operation performed: add Kg to the rLefsHandSideMatrix
  this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  // operation performed: add Kup to the rLefsHandSideMatrix
  this->CalculateAndAddKuJ( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  // operation performed: add Kpu to the rLefsHandSideMatrix
  this->CalculateAndAddKJu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  // operation performed: add Kpp to the rLefsHandSideMatrix
  this->CalculateAndAddKJJ( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  // operation performed: add Kpp Stab to the rLefsHandSideMatrix
  this->CalculateAndAddKJJStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  rVariables.detF = DeterminantF;
  rVariables.detF0 /= rVariables.detF;
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{

  rVariables.detF0 *= rVariables.detF;
  double DeterminantF = rVariables.detF;
  rVariables.detF = 1.0;

  //contribution of the internal and external forces
  VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

  // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
  this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

  // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
  this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );

  // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
  this->CalculateAndAddJacobianForces( rRightHandSideVector, rVariables, rIntegrationWeight );

  // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
  this->CalculateAndAddStabilizedJacobian( rRightHandSideVector, rVariables, rIntegrationWeight );

  rVariables.detF = DeterminantF;
  rVariables.detF0 /= rVariables.detF;
}

//************************** INTERNAL FORCES    *******************************
//************************************** Idem but with Total Stress ***********

void UpdatedLagrangianUJElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
                                                               ElementDataType & rVariables,
                                                               double& rIntegrationWeight)
{
  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().PointsNumber();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType node_dofs = this->GetNodeDofsSize();

  VectorType Fh=rRightHandSideVector;

  Vector StressVector = rVariables.StressVector;

  Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );

  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    SizeType indexup = node_dofs * i;
    SizeType indexu  = dimension * i;

    for ( SizeType j = 0; j < dimension; j++ )
    {
      rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
    }
  }

  KRATOS_CATCH( "" )
}

//******************************** JACOBIAN FORCES  **********************************
//************************************************************************************
void UpdatedLagrangianUJElement::CalculateAndAddJacobianForces(VectorType& rRightHandSideVector,
                                                               ElementDataType & rVariables,
                                                               double& rIntegrationWeight)
{
  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().PointsNumber();
  SizeType dimension = GetGeometry().WorkingSpaceDimension();

  SizeType indexp = dimension;

  VectorType Fh=rRightHandSideVector;


  double consistent;
  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    for ( SizeType j = 0; j < number_of_nodes; j++ )
    {
      if ( dimension == 2) {
        consistent = 1.0/12.0;
      }
      else {
        consistent = 1.0/20.0;
      }
      if ( i == j)
        consistent *= 2.0;

      const double& rNodalJacobian = (GetGeometry()[j].GetSolutionStepValue( JACOBIAN) );

      rRightHandSideVector[indexp] +=   consistent  * rNodalJacobian * rIntegrationWeight / Variables.detF0;

    }

    rRightHandSideVector[indexp] -= rVariables.N[i] * rIntegrationWeight;

    indexp += (dimension + 1);

  }

  KRATOS_CATCH( "" )

}

//******************************** STABILIZATION *************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddStabilizedJacobian(VectorType& rRightHandSideVector,
                                                                   ElementDataType & rVariables,
                                                                   double& rIntegrationWeight)
{
  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().PointsNumber();
  SizeType dimension = GetGeometry().WorkingSpaceDimension();

  SizeType indexp = dimension;

  //use of this variable for the complete parameter: (deffault: 4)
  double AlphaStabilization  = 4.0;

  //stabilization factor
  double StabilizationFactor = 0.20;
  if( GetProperties().Has(STABILIZATION_FACTOR_J) ){
    StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
  }
  else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR_J) ){
    StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR_J];
  }

  AlphaStabilization *= StabilizationFactor;

  const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
  const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

  double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
  double BulkModulus= YoungModulus/(3*(1-2*PoissonCoefficient));

  AlphaStabilization=(AlphaStabilization/(LameMu));

  AlphaStabilization *= BulkModulus;

  double ElementStabilizationFactor = 1.0;
  this->CalculateStabilizationFactor(ElementStabilizationFactor);

  if (YoungModulus < 0.00001)
  {
    AlphaStabilization = 4.0 * StabilizationFactor ;
    AlphaStabilization *= ElementStabilizationFactor;
  }

  if ( dimension == 2) {
    AlphaStabilization /= 18.0;
  }
  else {
    AlphaStabilization /= 80.0;
  }


  double consistent = 1;

  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    for ( SizeType j = 0; j < number_of_nodes; j++ )
    {

      consistent=(-1.0)*AlphaStabilization;
      if(i==j) {
        consistent=2.0*AlphaStabilization;
        if ( dimension == 3)
          consistent = 3.0 * AlphaStabilization;
      }

      double& Jacobian= GetGeometry()[j].FastGetSolutionStepValue(JACOBIAN);
      rRightHandSideVector[indexp] += consistent * Jacobian * rIntegrationWeight / rVariables.detF0;

    }

    indexp += (dimension + 1);
  }

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
//It includes the pw geometric stiffness

void UpdatedLagrangianUJElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                                     ElementDataType& rVariables,
                                                     double& rIntegrationWeight)
{
  KRATOS_TRY

  //assemble into rk the material uu contribution:
  const SizeType number_of_nodes = GetGeometry().PointsNumber();
  SizeType dimension = GetGeometry().WorkingSpaceDimension();
  double dimension_double = double(dimension);

  Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;

  SizeType voigtsize = 3;
  if (dimension == 3)
    voigtsize = 6;

  Matrix DeviatoricTensor(voigtsize,voigtsize);
  noalias( DeviatoricTensor) = ZeroMatrix(voigtsize,voigtsize);
  Vector Identity(voigtsize);
  noalias( Identity) = ZeroVector(voigtsize);
  for (SizeType i = 0; i < voigtsize ; ++i) {
    DeviatoricTensor(i,i) = 1.0;
  }
  for (SizeType i = 0; i < dimension; i++) {
    Identity(i) = 1.0;
    for (SizeType j = 0; j < dimension; j++) {
      DeviatoricTensor( i,j) -= 1.0/dimension_double;
    }
  }


  ConstitutiveMatrix = prod( ConstitutiveMatrix, DeviatoricTensor);

  Matrix AuxMatrix(voigtsize,voigtsize);
  noalias( AuxMatrix ) = ZeroMatrix(voigtsize,voigtsize);

  for (SizeType i = 0; i < voigtsize; i++) {
    for (SizeType j = 0; j < voigtsize; j++) {
      ConstitutiveMatrix(i,j)  += (1 -  2/dimension_double) * rVariables.StressVector(i) * Identity(j);
      AuxMatrix(i,j) += rVariables.StressVector(i) * Identity(j);
    }
  }

  const SizeType MatSize = dimension*number_of_nodes;
  MatrixType Kuu(MatSize,MatSize);

  noalias( Kuu ) = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( ConstitutiveMatrix, rVariables.B ) ) );


  SizeType indexi = 0;
  SizeType indexj  = 0;
  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    for ( SizeType idim = 0; idim < dimension ; idim ++)
    {
      indexj=0;
      for ( SizeType j = 0; j < number_of_nodes; j++ )
      {
        for ( SizeType jdim = 0; jdim < dimension ; jdim ++)
        {
          rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
          indexj++;
        }
      }
      indexi++;
    }
  }


  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddKuJ(MatrixType& rLeftHandSideMatrix,
                                                    ElementDataType& rVariables,
                                                    double& rIntegrationWeight)
{

  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  double dimension_double = double(dimension);

  Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;
  SizeType voigtsize = 3;
  if (dimension == 3)
    voigtsize = 6;


  // Trying to do it new
  Vector Identity(voigtsize);
  noalias( Identity) = ZeroVector(voigtsize);
  for (SizeType i = 0; i < dimension; i++)
    Identity(i) = 1.0;

  Vector ConstVector(voigtsize);
  noalias( ConstVector) = prod( ConstitutiveMatrix, Identity);
  ConstVector /= dimension_double;

  ConstVector += ( 2.0/dimension_double-1.0) * rVariables.StressVector;

  double ElementJacobian = 0.0;

  for ( SizeType i = 0; i <  number_of_nodes ; i++)
    ElementJacobian += GetGeometry()[i].GetSolutionStepValue( JACOBIAN ) * rVariables.N[i] ;

  ConstVector /= ElementJacobian;

  Vector KuJ(number_of_nodes*dimension);
  noalias( KuJ ) = prod( trans( rVariables.B), (ConstVector) );

  const SizeType MatSize = dimension*number_of_nodes;
  Matrix SecondMatrix(MatSize,MatSize);
  noalias(  SecondMatrix ) = ZeroMatrix( dimension*number_of_nodes, number_of_nodes);

  for (SizeType i = 0; i < dimension*number_of_nodes; i++) {
    for (SizeType j = 0; j < number_of_nodes; j++) {
      SecondMatrix(i,j) =KuJ(i) * rVariables.N[j];
    }
  }
  SecondMatrix *= rIntegrationWeight;


  // add the matrix in its place
  for (SizeType i = 0; i < number_of_nodes; i++) {
    for (SizeType idim = 0; idim < dimension; idim++) {
      for (SizeType j = 0; j < number_of_nodes; j++) {
        rLeftHandSideMatrix(i*(dimension+1) + idim, (dimension+1)*(j+1) -1 ) += SecondMatrix( i*(dimension) + idim, j);
      }
    }
  }

  KRATOS_CATCH( "" )
}



//******************* Kuug ********************************************************
//*********************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                                     ElementDataType& rVariables,
                                                     double& rIntegrationWeight)

{
  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  int size = number_of_nodes * dimension;


  Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

  Matrix ReducedKg = prod( rVariables.DN_DX,  rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized

  Matrix Kuu = zero_matrix<double> (size);
  MathUtils<double>::ExpandAndAddReducedMatrix( Kuu, ReducedKg, dimension );

  // MatrixType Kh=rLeftHandSideMatrix;

  //assemble into rLeftHandSideMatrix the geometric uu contribution:
  SizeType indexi = 0;
  SizeType indexj = 0;
  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    for ( SizeType idim = 0; idim < dimension ; idim ++)
    {
      indexj=0;
      for ( SizeType j = 0; j < number_of_nodes; j++ )
      {
        for ( SizeType jdim = 0; jdim < dimension ; jdim ++)
        {
          rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
          indexj++;
        }
      }
      indexi++;
    }
  }

  // std::cout<<std::endl;
  // std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddKJu(MatrixType& rLeftHandSideMatrix,
                                                    ElementDataType& rVariables,
                                                    double& rIntegrationWeight)

{
  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  MatrixType Kh=rLeftHandSideMatrix;

  //contributions to stiffness matrix calculated on the reference configuration
  SizeType indexp = dimension;


  for (SizeType i = 0; i < number_of_nodes; i++)
  {
    for (SizeType j = 0; j < number_of_nodes; j++)
    {
      int indexup = dimension*j + j;
      for (SizeType k = 0; k < dimension; k++)
      {
        rLeftHandSideMatrix(indexp, indexup + k) += rVariables.N[i] * rVariables.DN_DX( j, k) * rIntegrationWeight ;
      }
    }
    indexp += (dimension + 1);
  }

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddKJJ(MatrixType& rLeftHandSideMatrix,
                                                    ElementDataType& rVariables,
                                                    double& rIntegrationWeight)
{
  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  Matrix TotalF = prod( rVariables.F, rVariables.F0);

  MatrixType Kh=rLeftHandSideMatrix;

  //contributions to stiffness matrix calculated on the reference configuration
  SizeType indexpi = dimension;
  double consistent;

  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    SizeType indexpj = dimension;
    for ( SizeType j = 0; j < number_of_nodes; j++ )
    {
      if ( dimension == 2) {
        consistent = 1.0/12.0;
      }
      else {
        consistent = 1.0/20.0;
      }
      if ( i == j)
        consistent *= 2.0;

      rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * rIntegrationWeight / rVariables.detF0;
      indexpj += (dimension+1);
    }

    indexpi += (dimension + 1);
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateAndAddKJJStab(MatrixType& rLeftHandSideMatrix,
                                                        ElementDataType & rVariables,
                                                        double& rIntegrationWeight)
{

  KRATOS_TRY

  //repasar

  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  // MatrixType Kh=rLeftHandSideMatrix;

  //contributions to stiffness matrix calculated on the reference configuration
  SizeType indexpi = dimension;
  double consistent = 1.0;

  //use of this variable for the complete parameter: (deffault: 4)
  double AlphaStabilization  = 4.0;
  //stabilization factor
  double StabilizationFactor = 0.20;
  if( GetProperties().Has(STABILIZATION_FACTOR_J) ){
    StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
  }
  else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR_J) ){
    StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR_J];
  }

  AlphaStabilization *= StabilizationFactor;

  const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
  const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

  double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
  double BulkModulus= YoungModulus/(3*(1-2*PoissonCoefficient));

  AlphaStabilization=(AlphaStabilization/(LameMu));

  AlphaStabilization *= BulkModulus;  // TIMES THE BULK MODULUS BECAUSE I HAVE ALL THE EQUATION MULTIPLIED BY THE BULK MODULUS

  double ElementStabilizationFactor = 1.0;
  this->CalculateStabilizationFactor(ElementStabilizationFactor);

  if (YoungModulus < 0.00001)
  {
    AlphaStabilization = 4.0 * StabilizationFactor ;
    AlphaStabilization *= ElementStabilizationFactor;
  }

  if ( dimension == 2) {
    AlphaStabilization /= 18.0;
  }
  else {
    AlphaStabilization /= 80.0;
  }

  for ( SizeType i = 0; i < number_of_nodes; i++ )
  {
    SizeType indexpj = dimension;
    for ( SizeType j = 0; j < number_of_nodes; j++ )
    {
      consistent=(-1.0)*AlphaStabilization;
      if(i==j) {
        consistent=2.0*AlphaStabilization;
        if ( dimension == 3)
          consistent = 3.0 * AlphaStabilization;
      }

      rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

      indexpj += (dimension + 1);
    }

    indexpi += (dimension + 1);
  }

  // std::cout<<std::endl;
  // std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;
  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType MatSize = this->GetDofsSize();
    const SizeType node_dofs = this->GetNodeDofsSize();

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

    noalias(rMassMatrix) = ZeroMatrix( MatSize, MatSize );

    // Not Lumped Mass Matrix (numerical integration):

    //reading integration points
    IntegrationMethod CurrentIntegrationMethod = mThisIntegrationMethod; //GeometryData::GI_GAUSS_2; //GeometryData::GI_GAUSS_1;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( CurrentIntegrationMethod  );

    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);


    for ( SizeType PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
      //compute element kinematics
      this->CalculateKinematics( Variables, PointNumber );

      //getting informations for integration
      Variables.IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;

      Variables.IntegrationWeight = this->CalculateIntegrationWeight( Variables.IntegrationWeight );

      //compute point volume change
      double PointVolumeChange = 0;
      PointVolumeChange = this->CalculateVolumeChange( PointVolumeChange, Variables );

      double CurrentDensity = PointVolumeChange * GetProperties()[DENSITY];

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      	{
      	  SizeType indexupi = node_dofs * i;

      	  for ( SizeType j = 0; j < number_of_nodes; j++ )
      	    {
      	      SizeType indexupj = node_dofs * j;

      	      for ( SizeType k = 0; k < dimension; k++ )
      		{
      		  rMassMatrix( indexupi+k , indexupj+k ) += Variables.N[i] * Variables.N[j] * CurrentDensity * Variables.IntegrationWeight;
      		}
      	    }
      	}

    }

    // Lumped Mass Matrix:

    // double TotalMass = 0;

    // this->CalculateTotalMass( TotalMass, rCurrentProcessInfo );

    // if ( dimension == 2 ){
    //   if ( this->GetProperties().Has( THICKNESS ) )
    // 	TotalMass *= GetProperties()[THICKNESS];
    // }

    // Vector LumpFact(number_of_nodes);
    // noalias(LumpFact) = ZeroVector(number_of_nodes);

    // LumpFact = GetGeometry().LumpingFactors( LumpFact );

    // for ( SizeType i = 0; i < number_of_nodes; i++ )
    // {
    //     double temp = LumpFact[i] * TotalMass;

    //     SizeType indexup = node_dofs * i;

    //     for ( SizeType j = 0; j < dimension; j++ )
    //     {
    //         rMassMatrix( indexup+j , indexup+j ) = temp;
    //     }
    // }

    //std::cout<<std::endl;
    //std::cout<<" Mass Matrix "<<rMassMatrix<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  //0.-Initialize the DampingMatrix:
  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType MatSize = this->GetDofsSize();
  const SizeType node_dofs = this->GetNodeDofsSize();

  if ( rDampingMatrix.size1() != MatSize )
    rDampingMatrix.resize( MatSize, MatSize, false );

  noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );


  //1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
  double alpha = 0;
  if( GetProperties().Has(RAYLEIGH_ALPHA) ){
    alpha = GetProperties()[RAYLEIGH_ALPHA];
  }
  else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ){
    alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
  }

  double beta  = 0;
  if( GetProperties().Has(RAYLEIGH_BETA) ){
    beta = GetProperties()[RAYLEIGH_BETA];
  }
  else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ){
    beta = rCurrentProcessInfo[RAYLEIGH_BETA];
  }

  if( alpha != 0 || beta != 0){

    //1.-Calculate StiffnessMatrix:

    MatrixType LHSMatrix  = Matrix();

    this->CalculateLeftHandSide( LHSMatrix, rCurrentProcessInfo );

    MatrixType StiffnessMatrix  = Matrix();

    if ( StiffnessMatrix.size1() != MatSize )
      StiffnessMatrix.resize( MatSize, MatSize, false );

    StiffnessMatrix = ZeroMatrix( MatSize, MatSize );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      SizeType indexup = node_dofs * i;

      for ( SizeType j = 0; j < dimension; j++ )
      {
        StiffnessMatrix( indexup+j , indexup+j ) = LHSMatrix( indexup+j , indexup+j );
      }
    }

    //2.-Calculate MassMatrix:

    MatrixType MassMatrix  = Matrix();

    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );

    //3.-Compose the Damping Matrix:

    //Rayleigh Damping Matrix: alpha*M + beta*K
    rDampingMatrix  = alpha * MassMatrix;
    rDampingMatrix += beta  * StiffnessMatrix;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber )
{
  KRATOS_TRY

  LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

  //Deformation Gradient F0
  rVariables.detF0 = mDeterminantF0[rPointNumber];
  rVariables.F0    = mDeformationGradientF0[rPointNumber];

  rVariables.detH = mDeterminantJ0[rPointNumber];
  rVariables.H    = mDeformationGradientJ0[rPointNumber];


  KRATOS_CATCH("")
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianUJElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
{
  KRATOS_TRY

  rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

  return rVolumeChange;

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::SetElementData(ElementDataType& rVariables,
                                                ConstitutiveLaw::Parameters& rValues,
                                                const int & rPointNumber)
{
  KRATOS_TRY

  //LM is needed ?
  //rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

  //set previous step for output print purposes
  if( this->Is(SolidElement::FINALIZED_STEP) ){
    this->GetHistoricalVariables(rVariables, rPointNumber);
  }

  //check inverted element
  this->CheckElementData(rVariables, rPointNumber);

  rValues.SetStrainVector(rVariables.StrainVector);
  rValues.SetStressVector(rVariables.StressVector);
  rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
  rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
  rValues.SetShapeFunctionsValues(rVariables.N);

  ////
  //calculate nodal deformation gradient
  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = GetGeometry().size();

  // Compute F and detF (from 0 to n+1) : store it in H variable and detH
  rVariables.detH = rVariables.detF * rVariables.detF0;
  noalias(rVariables.H) = prod( rVariables.F, rVariables.F0 );

  // Calculate integration point JACOBIAN (total deformation gradient)
  double detJ = 0;
  for(SizeType i=0; i<number_of_nodes; ++i)
    detJ += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue(JACOBIAN);

  //add the effect of the interpolation
  double power = 1.0/double(dimension);
  rVariables.H *= pow(detJ/rVariables.detH, power);
  rVariables.detH = detJ;

  double detF0;
  Matrix invF0;
  MathUtils<double>::InvertMatrix(rVaribales.F0, invF0, detF0);

  if( !this->Is(SolidElement::FINALIZED_STEP) ){
    detJ = 0;
    for(SizeType i=0; i<number_of_nodes; ++i)
      detJ += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue(JACOBIAN,1);
  }

  invF0 /= pow(detJ/rVariables.detF0, power);

  Matrix F;
  noalias(F)= prod(rVariables.H, invF0);
  noalias(rVariables.H) = prod(F, mDeformationGradientJ0[rPointNumber]);
  //calculate nodal deformation gradient
  ////

  //set deformation gradient
  rValues.SetDeterminantF(rVariables.detJ);
  rValues.SetDeformationGradientF(rVariables.H);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

int UpdatedLagrangianUJElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  // Perform base element checks
  int ErrorCode = 0;
  ErrorCode = LargeDisplacementElement::Check(rCurrentProcessInfo);

  // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
  for(SizeType i=0; i<this->GetGeometry().size(); ++i)
  {
    // Nodal data
    Node<3> &rNode = this->GetGeometry()[i];
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rNode);
    //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);

    // Nodal dofs
    KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rNode);
    KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rNode);
    if( rCurrentProcessInfo[SPACE_DIMENSION] == 3)
      KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rNode);
  }

  // Check compatibility with the constitutive law
  ConstitutiveLaw::Features LawFeatures;
  this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

  if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
    KRATOS_ERROR << "constitutive law is not compatible with the U-J element type " << std::endl;

  // Check that all required variables have been registered
  KRATOS_CHECK_VARIABLE_KEY(JACOBIAN);

  return ErrorCode;

  KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUJElement::save( Serializer& rSerializer ) const
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
  rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
  rSerializer.save("DeterminantF0",mDeterminantF0);
  rSerializer.save("DeformationGradientJ0",mDeformationGradientJ0);
  rSerializer.save("DeterminantJ0",mDeterminantJ0);

}

void UpdatedLagrangianUJElement::load( Serializer& rSerializer )
{
  KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
  rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
  rSerializer.load("DeterminantF0",mDeterminantF0);
  rSerializer.load("DeformationGradientJ0",mDeformationGradientJ0);
  rSerializer.load("DeterminantJ0",mDeterminantJ0);
}

}
