// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//                   Manuel Caicedo
//                   Alfredo Huespe
//  Collaborator:    Vicente Mataix Ferrandiz
//

#include "linear_j2_plasticity_plane_strain_2d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::LinearJ2PlasticityPlaneStrain2D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::LinearJ2PlasticityPlaneStrain2D(const LinearJ2PlasticityPlaneStrain2D &rOther)
            : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearJ2PlasticityPlaneStrain2D::Clone() const
{
    LinearJ2PlasticityPlaneStrain2D::Pointer pclone(new LinearJ2PlasticityPlaneStrain2D(*this));
    return pclone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::~LinearJ2PlasticityPlaneStrain2D()
{
}

//************************************************************************************
//************************************************************************************

bool LinearJ2PlasticityPlaneStrain2D::Has(const Variable<bool>& rThisVariable)
{
    if(rThisVariable == INELASTIC_FLAG){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

bool LinearJ2PlasticityPlaneStrain2D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == PLASTIC_STRAIN){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

bool& LinearJ2PlasticityPlaneStrain2D::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

double& LinearJ2PlasticityPlaneStrain2D::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == PLASTIC_STRAIN){
        rValue = mAccumulatedPlasticStrain;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rThisVariable == INELASTIC_FLAG){
        mInelasticFlag = rValue;
    }
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rThisVariable == PLASTIC_STRAIN){
        mAccumulatedPlasticStrain = rValue;
    }
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::InitializeMaterial(const Properties& material_prop,
                                                         const GeometryType& rElementGeometry,
                                                         const Vector& rShapeFunctionsValues)
{
    mPlasticStrainOld = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrainOld = 0.0;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    mPlasticStrainOld = mPlasticStrain;
    mAccumulatedPlasticStrainOld = mAccumulatedPlasticStrain;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    Flags &Options = rValues.GetOptions();

    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    Vector& strain_vector = rValues.GetStrainVector();
    Vector& stress_vector = rValues.GetStressVector();

    if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(ConstitutiveMatrix, rMaterialProperties);
    }

    if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        Matrix elastic_tensor;
        Matrix& tangent_tensor = rValues.GetConstitutiveMatrix();
        const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
        const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
        const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
        const double E = rMaterialProperties[YOUNG_MODULUS];
        const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
        const double mu = E / (2. + 2. * poisson_ratio);
        const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));
        const double sqrt_two_thirds = std::sqrt(2.0 / 3.0); // =0.8164965809277260
        double trial_yield_function;

        mPlasticStrain = mPlasticStrainOld;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;

        elastic_tensor.resize(4, 4, false);
        CalculateElasticMatrix(elastic_tensor, rMaterialProperties);
        Vector sigma_trial(4);
        noalias(sigma_trial) = prod(elastic_tensor, strain_vector - mPlasticStrainOld);

        // stress_trial_dev = sigma - 1/3 tr(sigma) * I
        Vector stress_trial_dev = sigma_trial;

        const double trace = 1.0 / 3.0 * (sigma_trial(0) + sigma_trial(1) + sigma_trial(2));
        stress_trial_dev(0) -= trace;
        stress_trial_dev(1) -= trace;
        stress_trial_dev(2) -= trace;
        const double norm_dev_stress = std::sqrt(stress_trial_dev(0) * stress_trial_dev(0) +
                                           stress_trial_dev(1) * stress_trial_dev(1) +
                                           stress_trial_dev(2) * stress_trial_dev(2) +
                                           2. * stress_trial_dev(3) * stress_trial_dev(3));
        trial_yield_function = this->YieldFunction(norm_dev_stress, rMaterialProperties);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            mInelasticFlag = false;
            stress_vector = sigma_trial;
            tangent_tensor = elastic_tensor;
        }
        else {
            // INELASTIC
            mInelasticFlag = true;
            double dgamma = 0;
            Vector yield_function_normal_vector = stress_trial_dev / norm_dev_stress;
            if (delta_k != 0.0 && hardening_exponent != 0.0) {
                // Exponential softening
                dgamma = GetDeltaGamma(norm_dev_stress, rMaterialProperties);
            }
            else {
                // Linear softening
                dgamma = trial_yield_function /
                         (2. * mu * (1. + (hardening_modulus / (3. * mu))));
            }

            stress_vector(0) =
                volumetric_modulus * (strain_vector(0) + strain_vector(1) + strain_vector(2)) +
                stress_trial_dev(0) - 2. * mu * dgamma * yield_function_normal_vector(0);
            stress_vector(1) =
                volumetric_modulus * (strain_vector(0) + strain_vector(1) + strain_vector(2)) +
                stress_trial_dev(1) - 2. * mu * dgamma * yield_function_normal_vector(1);
            stress_vector(2) =
                volumetric_modulus * (strain_vector(0) + strain_vector(1) + strain_vector(2)) +
                stress_trial_dev(2) - 2. * mu * dgamma * yield_function_normal_vector(2);
            stress_vector(3) =
                stress_trial_dev(3) - 2. * mu * dgamma * yield_function_normal_vector(3);

            mPlasticStrain(0) = mPlasticStrainOld(0) + dgamma * yield_function_normal_vector(0);
            mPlasticStrain(1) = mPlasticStrainOld(1) + dgamma * yield_function_normal_vector(1);
            mPlasticStrain(2) = mPlasticStrainOld(2) + dgamma * yield_function_normal_vector(2);
            mPlasticStrain(3) = mPlasticStrainOld(3) + dgamma * yield_function_normal_vector(3) * 2;
            mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + sqrt_two_thirds * dgamma;

            // Update derivative of the hardening-softening modulus

            CalculateTangentTensor(dgamma, norm_dev_stress, yield_function_normal_vector,
                                   rMaterialProperties, tangent_tensor);
        }
    }
}

//************************************************************************************
//************************************************************************************

double& LinearJ2PlasticityPlaneStrain2D::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == STRAIN_ENERGY){
        Flags &Options = rParameterValues.GetOptions();
        Options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        CalculateMaterialResponseCauchy(rParameterValues);
        Vector& strain_vector = rParameterValues.GetStrainVector();
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        Matrix elastic_tensor(6, 6);
        CalculateElasticMatrix(elastic_tensor, r_material_properties);
        // Linear + exponential hardening
        rValue = 0.5 * inner_prod(strain_vector - mPlasticStrain, prod(elastic_tensor, strain_vector - mPlasticStrain)) + GetPlasticPotential(r_material_properties);
    }
    if(rThisVariable == PLASTIC_STRAIN){
        rValue = mAccumulatedPlasticStrain;
    }
    return(rValue);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponsePK1(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************


void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************


void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

double LinearJ2PlasticityPlaneStrain2D::GetSaturationHardening(const Properties& rMaterialProperties)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

    const double k_new = yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrain) +
               delta_k * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
    return k_new;
}

//************************************************************************************
//************************************************************************************

double LinearJ2PlasticityPlaneStrain2D::GetPlasticPotential(const Properties& rMaterialProperties)
{
   const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
   const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
   const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
   const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

   const double wp_new = 0.5*(theta * hardening_modulus * std::pow(mAccumulatedPlasticStrain, 2.0)) +
                   delta_k * (mAccumulatedPlasticStrain -
                   (1/hardening_exponent) * (1- std::exp(-hardening_exponent * mAccumulatedPlasticStrain)));
   return wp_new;
}

//************************************************************************************
//************************************************************************************

double LinearJ2PlasticityPlaneStrain2D::GetDeltaGamma(
    const double NormStressTrial,
    const Properties& rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double tolerance = 1e-6 * yield_stress;
    const double mu = E / (2. * (1. + poisson_ratio));
    const double sqrt_two_thirds = std::sqrt(2.0 / 3.0); // =0.8164965809277260
    double dgamma = 0.0;
    double norm_yieldfunction = 1.0;

    while (norm_yieldfunction > tolerance)
    {
        const double k_new = GetSaturationHardening(rMaterialProperties);
        const double kp_new = theta * hardening_modulus +
            delta_k * (hardening_exponent * std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
        const double yieldfunction = - sqrt_two_thirds * k_new + NormStressTrial - 2. * mu * dgamma;
        const double derivative_yieldfunction = -2. * mu * (1. + kp_new / (3. * mu));
        dgamma = dgamma - yieldfunction / derivative_yieldfunction;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + sqrt_two_thirds * dgamma;
        norm_yieldfunction = std::abs(yieldfunction);
    }
    // TODO (marcelo): handle the case when no convergence is achieved.
    return dgamma;
}

//************************************************************************************
//************************************************************************************

double LinearJ2PlasticityPlaneStrain2D::YieldFunction(
    const double NormDeviationStress,
    const Properties& rMaterialProperties
    )
{
    const double sqrt_two_thirds = std::sqrt(2.0 / 3.0);
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double k_old =
        yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrainOld) +
        (delta_k) * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrainOld));

    return NormDeviationStress - k_old * sqrt_two_thirds;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateElasticMatrix(
    Matrix &rElasticityTensor,
    const Properties &rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (rElasticityTensor.size1() != 4 || rElasticityTensor.size2() != 4)
        rElasticityTensor.resize(4, 4, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2. * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(0, 3) = 0.;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2. * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(1, 3) = 0.;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2. * mu;
    rElasticityTensor(2, 3) = 0.;
    rElasticityTensor(3, 0) = 0.;
    rElasticityTensor(3, 1) = 0.;
    rElasticityTensor(3, 2) = 0.;
    rElasticityTensor(3, 3) = mu;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateTangentTensor(
    const double DeltaGamma,
    const double NormStressTrial,
    const Vector& YieldFunctionNormalVector,
    const Properties& rMaterialProperties,
    Matrix& rElasticityTensor
    )
{
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double mu = E / (2. + 2. * poisson_ratio);
    const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

    const double kp_new = (theta * hardening_modulus) +
                    delta_k * (hardening_exponent *
                               std::exp(-hardening_exponent * mAccumulatedPlasticStrain));

    const double theta_new = 1 - (2. * mu * DeltaGamma) / NormStressTrial;
    const double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

    rElasticityTensor(0, 0) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(0)));
    rElasticityTensor(0, 1) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(1)));
    rElasticityTensor(0, 2) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(2)));
    rElasticityTensor(0, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(3)));

    rElasticityTensor(1, 0) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(0)));
    rElasticityTensor(1, 1) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(1)));
    rElasticityTensor(1, 2) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(2)));
    rElasticityTensor(1, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(3)));

    rElasticityTensor(2, 0) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(0)));
    rElasticityTensor(2, 1) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(1)));
    rElasticityTensor(2, 2) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(2)));
    rElasticityTensor(2, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(3)));

    rElasticityTensor(3, 0) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(0)));
    rElasticityTensor(3, 1) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(1)));
    rElasticityTensor(3, 2) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(2)));
    rElasticityTensor(3, 3) = mu * theta_new - (2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(3)));
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 4;
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

int LinearJ2PlasticityPlaneStrain2D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(REFERENCE_HARDENING_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(INFINITY_HARDENING_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_EXPONENT));

    return 0;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mPlasticStrainOld", mPlasticStrainOld);
    rSerializer.save("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.save("mAccumulatedPlasticStrainOld", mAccumulatedPlasticStrainOld);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mInelasticFlag", mInelasticFlag);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mPlasticStrainOld", mPlasticStrainOld);
    rSerializer.load("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.load("mAccumulatedPlasticStrainOld", mAccumulatedPlasticStrainOld);
}

} /* namespace Kratos.*/
