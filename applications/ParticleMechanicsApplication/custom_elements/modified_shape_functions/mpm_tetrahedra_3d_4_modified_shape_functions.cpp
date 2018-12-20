//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//

// System includes

// External includes

// Project includes
#include "custom_elements/modified_shape_functions/mpm_tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{

/// MPMTetrahedra3D4ModifiedShapeFunctions implementation
/// Default constructor
MPMTetrahedra3D4ModifiedShapeFunctions::MPMTetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    MPMModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTetrahedraSplitter(Kratos::make_shared<DivideTetrahedra3D4>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTetrahedraSplitter->GenerateDivision();
    mpTetrahedraSplitter->GenerateIntersectionsSkin();
};

/// Destructor
MPMTetrahedra3D4ModifiedShapeFunctions::~MPMTetrahedra3D4ModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string MPMTetrahedra3D4ModifiedShapeFunctions::Info() const {
    return "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print information about this object.
void MPMTetrahedra3D4ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print object's data.
void MPMTetrahedra3D4ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Tetrahedra3D4N modified shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::stringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};

// Returns a pointer to the splitting utility
const DivideGeometry::Pointer MPMTetrahedra3D4ModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTetrahedraSplitter;
};

// Returns true if the element is splitting
bool MPMTetrahedra3D4ModifiedShapeFunctions::IsSplit() {
    return mpTetrahedraSplitter->mIsSplit;
};

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const array_1d<double,3>& rIntegrationPoint) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
                                    mpTetrahedraSplitter->mEdgeNodeI,
                                    mpTetrahedraSplitter->mEdgeNodeJ,
                                    mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side values
        this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues,
                                     rPositiveSideShapeFunctionsGradientsValues,
                                     rPositiveSideWeightsValues,
                                     mpTetrahedraSplitter->mPositiveSubdivisions,
                                     p_matrix,
                                     rIntegrationPoint);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative side.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const array_1d<double,3>& rIntegrationPoint) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
                                    mpTetrahedraSplitter->mEdgeNodeI,
                                    mpTetrahedraSplitter->mEdgeNodeJ,
                                    mpTetrahedraSplitter->mSplitEdges);

        // Compute the negative side values
        this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues,
                                     rNegativeSideShapeFunctionsGradientsValues,
                                     rNegativeSideWeightsValues,
                                     mpTetrahedraSplitter->mNegativeSubdivisions,
                                     p_matrix,
                                     rIntegrationPoint);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const array_1d<double,3>& rIntegrationPoint) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
                                    mpTetrahedraSplitter->mEdgeNodeI,
                                    mpTetrahedraSplitter->mEdgeNodeJ,
                                    mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rInterfacePositiveSideShapeFunctionsValues,
                                              rInterfacePositiveSideShapeFunctionsGradientsValues,
                                              rInterfacePositiveSideWeightsValues,
                                              mpTetrahedraSplitter->mPositiveInterfaces,
                                              mpTetrahedraSplitter->mPositiveSubdivisions,
                                              mpTetrahedraSplitter->mPositiveInterfacesParentIds,
                                              p_matrix,
                                              rIntegrationPoint);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const array_1d<double,3>& rIntegrationPoint) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
                                    mpTetrahedraSplitter->mEdgeNodeI,
                                    mpTetrahedraSplitter->mEdgeNodeJ,
                                    mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rInterfaceNegativeSideShapeFunctionsValues,
                                              rInterfaceNegativeSideShapeFunctionsGradientsValues,
                                              rInterfaceNegativeSideWeightsValues,
                                              mpTetrahedraSplitter->mNegativeInterfaces,
                                              mpTetrahedraSplitter->mNegativeSubdivisions,
                                              mpTetrahedraSplitter->mNegativeInterfacesParentIds,
                                              p_matrix,
                                              rIntegrationPoint);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Given a face id, computes the positive side subdivision shape function values in that face.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rPositiveExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
    Vector &rPositiveExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const array_1d<double,3>& rIntegrationPoint
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(
            p_matrix,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);

        // Get the external faces
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mPositiveSubdivisions,
            FaceId);

        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rPositiveExteriorFaceShapeFunctionsValues,
            rPositiveExteriorFaceShapeFunctionsGradientsValues,
            rPositiveExteriorFaceWeightsValues,
            exterior_faces_vector,
            mpTetrahedraSplitter->mPositiveSubdivisions,
            exterior_faces_parent_ids_vector,
            p_matrix,
            rIntegrationPoint);

    } else {
        KRATOS_ERROR << "Using the ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Given a face id, computes the positive side subdivision shape function values in that face.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rNegativeExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
    Vector &rNegativeExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const array_1d<double,3>& rIntegrationPoint
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(
            p_matrix,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);

        // Get the external faces
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mNegativeSubdivisions,
            FaceId);

        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rNegativeExteriorFaceShapeFunctionsValues,
            rNegativeExteriorFaceShapeFunctionsGradientsValues,
            rNegativeExteriorFaceWeightsValues,
            exterior_faces_vector,
            mpTetrahedraSplitter->mNegativeSubdivisions,
            exterior_faces_parent_ids_vector,
            p_matrix,
            rIntegrationPoint);

    } else {
        KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    std::vector<Vector> &rPositiveSideInterfaceAreaNormals,
    const array_1d<double,3>& rIntegrationPoint) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeFaceNormalOnOneSide(rPositiveSideInterfaceAreaNormals,
                                              mpTetrahedraSplitter->mPositiveInterfaces,
                                              rIntegrationPoint);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
    std::vector<Vector> &rNegativeSideInterfaceAreaNormals,
    const array_1d<double,3>& rIntegrationPoint) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeFaceNormalOnOneSide(rNegativeSideInterfaceAreaNormals,
                                              mpTetrahedraSplitter->mNegativeInterfaces,
                                              rIntegrationPoint);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// For a given face, computes the positive side exteriorface outwards area normal vector values.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputePositiveExteriorFaceAreaNormals(
    std::vector<Vector> &rPositiveExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const array_1d<double,3>& rIntegrationPoint)
{

    if (this->IsSplit())
    {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mPositiveSubdivisions,
            FaceId);

        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rPositiveExteriorFaceAreaNormal,
                                              exterior_faces_vector,
                                              rIntegrationPoint);
    }
    else
    {
        KRATOS_ERROR << "Using the ComputePositiveExteriorFaceAreaNormals method for a non divided geometry.";
    }
};

// For a given face, computes the positive side exterior face outwards area normal vector values.
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeExteriorFaceAreaNormals(
    std::vector<Vector> &rNegativeExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const array_1d<double,3>& rIntegrationPoint)
{

    if (this->IsSplit())
    {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mNegativeSubdivisions,
            FaceId);

        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rNegativeExteriorFaceAreaNormal,
                                              exterior_faces_vector,
                                              rIntegrationPoint);
    }
    else
    {
        KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceAreaNormals method for a non divided geometry.";
    }
};

// Computes the positive side shape function values in the edges intersections
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeShapeFunctionsOnPositiveEdgeIntersections(
    Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues){

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(
            p_matrix,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);

        // Compute the edge intersections shape function values
        this->ComputeEdgeIntersectionValuesOnOneSide(
            p_matrix,
            rPositiveEdgeIntersectionsShapeFunctionsValues);

    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnPositiveEdgeIntersections method for a non divided geometry.";
    }
};

// Computes the negative side shape function values in the edges intersections
void MPMTetrahedra3D4ModifiedShapeFunctions::ComputeShapeFunctionsOnNegativeEdgeIntersections(
    Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues){

    if (this->IsSplit()) {
        // Note that positive and negative sides values are equal for standard shape functions
        this->ComputeShapeFunctionsOnPositiveEdgeIntersections(
            rNegativeEdgeIntersectionsShapeFunctionsValues);
    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnNegativeEdgeIntersections method for a non divided geometry.";
    }
};

}; //namespace Kratos
