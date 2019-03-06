//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on Armin Geiser and Martin Fusseder work
//

#include "incompressible_adjoint_potential_flow_element.h"

namespace Kratos
{
    template <class TPrimalElement<int Dim, int NumNodes>>
    Element::Pointer IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleAdjointPotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    Element::Pointer IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleAdjointPotentialFlowElement(NewId, pGeom, pProperties));
        KRATOS_CATCH("");
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    Element::Pointer IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const 
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleAdjointPotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    Element::IntegrationMethod IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetIntegrationMethod() const 
    {
        return mPrimalElement.GetIntegrationMethod();
    }
  
    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::Initialize() 
    {   
        mPrimalElement.Initialize();
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
    {
        mPrimalElement.Data() = this->Data();
        mPrimalElement.Set(Flags(*this));
        mPrimalElement.InitializeSolutionStep(rCurrentProcessInfo);
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
    {

    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) 
    {
        mPrimalElement.CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) 
    {
        MatrixType tmp;
        mPrimalElement.CalculateLeftHandSide(tmp, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(tmp);                                    
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) 
    {
        mPrimalElement.CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) 
    {
        mPrimalElement.GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
    
    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) 
    {
        mPrimalElement.GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetValuesVector(Vector& rValues, int Step) 
    {
        KRATOS_TRY
        const IncompressibleAdjointPotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if (wake == 1) // wake element
        {
            if(rValues.size() != 2*TPrimalElement<int Dim, int NumNodes>::NumNodes)
                rValues.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);
            GetValuesOnSplitElement(rValues,distances);
            
        }else{ // normal element
            if(rValues.size() != NumNodes)
                rValues.resize(NumNodes, false);

            if(kutta == 0){
                for(unsigned int i=0; i<NumNodes; i++)
                    rValues[i] =GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            }else{
                for(unsigned int i=0; i<NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY;
        const double delta = this->GetPerturbationSize();
        ProcessInfo process_info = rCurrentProcessInfo;

        Vector RHS;
        Vector RHS_perturbed;

        mPrimalElement.CalculateRightHandSide(RHS, process_info);

        if (rOutput.size1() != NumNodes)
            rOutput.resize(Dim*NumNodes, RHS.size(), false);
        bool kutta_node_spotted=false;
        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            if(GetGeometry()[i_node].GetValue(TRAILING_EDGE)){
                kutta_node_spotted=true;
                break;
            }
        }

        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            for(unsigned int i_dim = 0; i_dim<Dim; i_dim++){
                if (GetGeometry()[i_node].Is(SOLID)){
                    mPrimalElement.GetGeometry()[i_node].GetInitialPosition()[i_dim] += delta;
                    mPrimalElement.GetGeometry()[i_node].Coordinates()[i_dim] += delta;

                    // compute LHS after perturbation
                    mPrimalElement.CalculateRightHandSide(RHS_perturbed, process_info);

                    //compute derivative of RHS w.r.t. design variable with finite differences
                    for(unsigned int i = 0; i < RHS.size(); ++i)
                        rOutput( (i_dim + i_node*Dim), i) = (RHS_perturbed[i] - RHS[i]) / delta;

                    // unperturb the design variable
                    mPrimalElement.GetGeometry()[i_node].GetInitialPosition()[i_dim] -= delta;
                    mPrimalElement.GetGeometry()[i_node].Coordinates()[i_dim] -= delta;
                }else{
                    for(unsigned int i = 0; i < RHS.size(); ++i)
                        rOutput( (i_dim + i_node*Dim), i) = 0.0;
                }
            }
        }
        if (kutta_node_spotted){ //remove kutta node lines
            for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
                if(GetGeometry()[i_node].GetValue(TRAILING_EDGE)){
                    for(unsigned int i_dim = 0; i_dim<Dim; i_dim++){
                        for(unsigned int i = 0; i < RHS.size(); ++i)
                                rOutput( (i_dim + i_node*Dim), i) = 0.0;
                    }
                }
            }  
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) 
    {
        const IncompressibleAdjointPotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if(wake == 0)//normal element
        {
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
                        rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                    else
                        rResult[i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
                }
            }
        }
        else//wake element
        {
            if (rResult.size() != 2*NumNodes)
                rResult.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }
        }


    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) 
    {
        const IncompressibleAdjointPotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if(wake == 0) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
                        rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        else//wake element
        {
            if (rElementalDofList.size() != 2*NumNodes)
                rElementalDofList.resize(2*NumNodes);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }
        }
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    int IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::Check(const ProcessInfo& rCurrentProcessInfo) 
    {

        KRATOS_TRY

        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "IncompressibleAdjointPotentialFlowElement found with Id 0 or negative", "")
        }

        if (this->GetGeometry().Area() <= 0)
        {
            std::cout << "error on IncompressibleAdjointPotentialFlowElement -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0", "")
        }

        for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
        {
            if (this->GetGeometry()[i].SolutionStepsDataHas(ADJOINT_VELOCITY_POTENTIAL) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing variable ADJOINT_VELOCITY_POTENTIAL on node ", this->GetGeometry()[i].Id())
        }

        return 0;

        KRATOS_CATCH("");
    }


    /// Turn back information as a string.
    template <class TPrimalElement<int Dim, int NumNodes>>
    std::string IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::Info() const 
    {
        std::stringstream buffer;
        buffer << "IncompressibleAdjointPotentialFlowElement #" << Id();
        return buffer.str();
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::PrintInfo(std::ostream& rOStream) const 
    {
        rOStream << "IncompressibleAdjointPotentialFlowElement #" << Id();
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::PrintData(std::ostream& rOStream) const 
    {
        pGetGeometry()->PrintData(rOStream);
    }

    /*PROTECTED*/

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetWakeDistances(array_1d<double,NumNodes>& distances)
    {
        noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
    {

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] > 0)
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }

        //negative part - sign is opposite to the previous case
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] < 0)
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
    }

    /*PRIVATE*/

    template <class TPrimalElement<int Dim, int NumNodes>>
    double IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::GetPerturbationSize()
    {
        const double delta = this->GetValue(SCALE_FACTOR);
        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
        return delta;
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::save(Serializer& rSerializer) const 
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    template <class TPrimalElement<int Dim, int NumNodes>>
    void IncompressibleAdjointPotentialFlowElement<TPrimalElement<int Dim, int NumNodes>>::load(Serializer& rSerializer) 
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Template class instantiation

    template class IncompressibleAdjointPotentialFlowElement<2, 3>;
} // namespace Kratos.

