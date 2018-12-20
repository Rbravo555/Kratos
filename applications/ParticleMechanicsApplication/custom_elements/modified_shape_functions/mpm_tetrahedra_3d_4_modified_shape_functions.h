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

#if !defined(KRATOS_MPM_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS)
#define      KRATOS_MPM_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "utilities/divide_tetrahedra_3d_4.h"
#include "custom_elements/modified_shape_functions/mpm_modified_shape_functions.h"

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

class MPMTetrahedra3D4ModifiedShapeFunctions : public MPMModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPMTetrahedra3D4ModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(MPMTetrahedra3D4ModifiedShapeFunctions);

    // General type definitions
    typedef MPMModifiedShapeFunctions                          BaseType;
    typedef BaseType::GeometryType                             GeometryType;
    typedef BaseType::GeometryPointerType                      GeometryPointerType;
    typedef BaseType::ShapeFunctionsGradientsType              ShapeFunctionsGradientsType;

    typedef BaseType::IndexedPointGeometryType                 IndexedPointGeometryType;
    typedef BaseType::IndexedPointGeometryPointerType          IndexedPointGeometryPointerType;

    typedef BaseType::IntegrationPointType                     IntegrationPointType;
    typedef BaseType::IntegrationPointsArrayType               IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    MPMTetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType rpInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    ~MPMTetrahedra3D4ModifiedShapeFunctions() override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * Returns the member pointer to the splitting utility.
    */
    const DivideGeometry::Pointer pGetSplittingUtil() const override;

    /**
    * Returns the shape function values in the positive split element side for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rPositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
        Vector &rPositiveSideWeightsValues,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputeNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
        Vector &rNegativeSideWeightsValues,
        const array_1d<double,3>& rIntegrationPoint) override;

    ///@}

    /**
    * Returns the shape function values in the positive split element interface side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfacePositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
        Vector &rInterfacePositiveSideWeightsValues,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Returns the shape function values in the negative split element interface side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfaceNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
        Vector &rInterfaceNegativeSideWeightsValues,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Given a face id, returns the shape function values in the positive split element exterior face side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param  FaceId Face local id. in where the values are to be computed.
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rPositiveExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
        Vector &rPositiveExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Given a face id, returns the shape function values in the negative split element exterior face side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param  FaceId Face local id. in where the values are to be computed.
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rNegativeExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
        Vector &rNegativeExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveSideInterfaceAreaNormals: Outwards area normal vector list.
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputePositiveSideInterfaceAreaNormals(
        std::vector<Vector> &rPositiveSideInterfaceAreaNormals,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeSideInterfaceAreaNormals: Outwards area normal vector list.
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputeNegativeSideInterfaceAreaNormals(
        std::vector<Vector> &rNegativeSideInterfaceAreaNormals,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param  FaceId Face local id. in where the values are to be computed.
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputePositiveExteriorFaceAreaNormals(
        std::vector<Vector> &rPositiveExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Returns the negative side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param  FaceId Face local id. in where the values are to be computed.
    * @param  rIntegrationPoint: Coordinate of integration point.
    */
    void ComputeNegativeExteriorFaceAreaNormals(
        std::vector<Vector> &rNegativeExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const array_1d<double,3>& rIntegrationPoint) override;

    /**
    * Returns the positive side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the positive side edge intersection shape function values. For non-split edges,
    * the corresponding row is plenty of zeros.
    */
    void ComputeShapeFunctionsOnPositiveEdgeIntersections(
        Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues) override;

    /**
    * Returns the negative side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the negative side edge intersection shape function values. For non-split edges,
    * the corresponding row is plenty of zeros.
    */
    void ComputeShapeFunctionsOnNegativeEdgeIntersections(
        Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues) override;

    /**
    * Returns true if the element is split and false otherwise.
    */
    bool IsSplit() override;

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

    DivideTetrahedra3D4::Pointer mpTetrahedraSplitter;

    ///@}
    ///@name Serialization
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

    /// Assignment operator.
    MPMTetrahedra3D4ModifiedShapeFunctions& operator=(MPMTetrahedra3D4ModifiedShapeFunctions const& rOther);

    /// Copy constructor.
    MPMTetrahedra3D4ModifiedShapeFunctions(MPMTetrahedra3D4ModifiedShapeFunctions const& rOther) :
        MPMModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()),
        mpTetrahedraSplitter(new DivideTetrahedra3D4(*rOther.GetInputGeometry(), rOther.GetNodalDistances())) {

        // Perform the element splitting
        mpTetrahedraSplitter->GenerateDivision();
        mpTetrahedraSplitter->GenerateIntersectionsSkin();
    };

    ///@}

};// class MPMTetrahedra3D4ModifiedShapeFunctions

}
#endif /* KRATOS_MPM_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS defined */
