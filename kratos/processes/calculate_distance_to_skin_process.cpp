//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


// System includes
#include <chrono>

// External includes


// Project includes
#include "processes/calculate_distance_to_skin_process.h"


namespace Kratos
{
namespace Internals
{
    class Timer {
    std::chrono::high_resolution_clock::time_point _start;

   public:
    using duration_type = long long;
    Timer() : _start(std::chrono::high_resolution_clock::now()) {}
    std::chrono::milliseconds elapsed() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - _start);
    }
    void reset() { _start = std::chrono::high_resolution_clock::now(); }
};
}
CalculateDistanceToSkinProcess::CalculateDistanceToSkinProcess(
    ModelPart &rVolumePart, ModelPart &rSkinPart, Parameters TheParameters)
    : CalculateDiscontinuousDistanceToSkinProcess(rVolumePart, rSkinPart) , mStatistics(){

  Parameters default_parameters(R"(
            {
                "epsilon"                       : 1e-12, 
                "ray_casting" : {
                    "xyz_valid_weights"       : [1.0, 1.0, 1.0],
                    "xyz_invalid_weights"     : [0.0, 0.0, 0.0]
                }
            }  )");

  TheParameters.ValidateAndAssignDefaults(default_parameters);

  mEpsilon = TheParameters["epsilon"].GetDouble();

  Parameters ray_casting_parameters = TheParameters["ray_casting"];
  mValidRayWeights[0] = ray_casting_parameters["xyz_valid_weights"][0].GetDouble();
  mValidRayWeights[1] = ray_casting_parameters["xyz_valid_weights"][1].GetDouble();
  mValidRayWeights[2] = ray_casting_parameters["xyz_valid_weights"][2].GetDouble();

  mInvalidRayWeights[0] = ray_casting_parameters["xyz_invalid_weights"][0].GetDouble();
  mInvalidRayWeights[1] = ray_casting_parameters["xyz_invalid_weights"][1].GetDouble();
  mInvalidRayWeights[2] = ray_casting_parameters["xyz_invalid_weights"][2].GetDouble();
}

CalculateDistanceToSkinProcess::~CalculateDistanceToSkinProcess() {}

void CalculateDistanceToSkinProcess::Initialize() {
  CalculateDiscontinuousDistanceToSkinProcess::Initialize();
  this->InitializeNodalDistances();
	}

	void CalculateDistanceToSkinProcess::InitializeNodalDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess).GetModelPart1();

		for (auto& node : ModelPart1.Nodes())
		{
			node.GetSolutionStepValue(DISTANCE) = std::numeric_limits<double>::max();
		}
	}

	void CalculateDistanceToSkinProcess::CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		CalculateDiscontinuousDistanceToSkinProcess::CalculateDistances(rIntersectedObjects);
		this->CalculateNodalDistances();
		this->CalculateNodesDistances();
	}

	void CalculateDistanceToSkinProcess::CalculateNodalDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess).GetModelPart1();

		constexpr int number_of_tetrahedra_points = 4;
		for (auto& element : ModelPart1.Elements())
		{
			if (element.Is(TO_SPLIT))
			{
				const auto& r_elemental_distances = element.GetValue(ELEMENTAL_DISTANCES);
				for (int i = 0; i < number_of_tetrahedra_points; i++)
				{
					Node<3>& r_node = element.GetGeometry()[i];
					double& r_distance = r_node.GetSolutionStepValue(DISTANCE);
					if (fabs(r_distance) > fabs(r_elemental_distances[i]))
						r_distance = r_elemental_distances[i];
				}
			}
		}
	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	void CalculateDistanceToSkinProcess::CalculateNodesDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess).GetModelPart1();

		#pragma omp parallel for
		for(int k = 0 ; k < static_cast<int>(ModelPart1.NumberOfNodes()); ++k)
        {
			ModelPart::NodesContainerType::iterator itNode = ModelPart1.NodesBegin() + k;
            this->CalculateNodeDistance(*itNode);
        }
        		
        constexpr int number_of_tetrahedra_points = 4;
        for(int repeat = 0 ; repeat < 10 ; repeat++){
            mStatistics.NumberOfIrregularElementalDistances = 0;
            for (auto& element : ModelPart1.Elements())
            {
                int number_of_positive_distances = 0;
                if (element.IsNot(TO_SPLIT))
                {
                    array_1d<double*, number_of_tetrahedra_points> distances;
                    for (int i = 0; i < number_of_tetrahedra_points; i++)
                    {
                        Node<3>& r_node = element.GetGeometry()[i];
                        distances[i] = &(r_node.GetSolutionStepValue(DISTANCE));
                        if(*(distances[i]) > -mEpsilon) 
                            number_of_positive_distances++;

                    }
                    if(number_of_positive_distances != 0 && number_of_positive_distances != number_of_tetrahedra_points){
                        if(number_of_positive_distances == 1){
                            for (int i = 0; i < number_of_tetrahedra_points; i++)
                                *(distances[i]) = -std::abs(*(distances[i]));
                        }
                        else if(number_of_positive_distances == 3){
                            for (int i = 0; i < number_of_tetrahedra_points; i++)
                                *(distances[i]) = std::abs(*(distances[i]));
                        }
                        else
                            mStatistics.NumberOfIrregularElementalDistances++;
                    }

                }
            }
            if(mStatistics.NumberOfIrregularElementalDistances == 0)
                break;
        }
	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	void CalculateDistanceToSkinProcess::CalculateNodeDistance(Node<3>& rNode)
	{
		double coord[3] = {rNode.X(), rNode.Y(), rNode.Z()};
		double distance = DistancePositionInSpace(coord);
		double& node_distance =  rNode.GetSolutionStepValue(DISTANCE);

		//const double epsilon = 1.00e-12;
		//if(fabs(node_distance) > fabs(distance))
		//    node_distance = distance;
		/*else*/ if (distance*node_distance < 0.00) // assigning the correct sign
			node_distance = -node_distance;
	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	double CalculateDistanceToSkinProcess::DistancePositionInSpace(double* coords)
	{

		typedef Element::GeometryType triangle_type;
        typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;

        intersections_container_type intersections;

        const int dimension = 3;

        double distances[3] = {1.0, 1.0, 1.0};
		bool ray_is_valid[3]={true,true,true};

        for (int i_direction = 0; i_direction < dimension; i_direction++)
        {
            // Creating the ray
            double ray[3] = {coords[0], coords[1], coords[2]};


			OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess::CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess.GetOctreePointer();
            pOctree->NormalizeCoordinates(ray);
            ray[i_direction] = 0; // starting from the lower extreme

            this->GetRayIntersections(ray, i_direction, intersections);

            int ray_color= 1;
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end())
			{
                double d = coords[i_direction] - i_intersection->first;
                if (d > mEpsilon)
				{
                    ray_color = -ray_color;
                    distances[i_direction] = d;
                }
				else if (d > -mEpsilon)
				{
                    distances[i_direction] = 0.00;
                    break;
                }
				else
				{
                    if(distances[i_direction] > -d)
                        distances[i_direction] = -d;
                    break;
                }
				ray_is_valid[i_direction] = !ray_is_valid[i_direction];
                i_intersection++;
            }
            if (ray_is_valid[i_direction])
                mStatistics.NumberOfValidRays[i_direction]++;
            else
                mStatistics.NumberOfInvalidRays[i_direction]++;

            distances[i_direction] *= ray_color;
        }

        double distance = (fabs(distances[0]) > fabs(distances[1])) ? distances[1] : distances[0];
        distance = (fabs(distance) > fabs(distances[2])) ? distances[2] : distance;

		// we vote
		double positive_distances_weight = 0;
		double negative_distances_weight = 0;
        bool has_positive_distance = false;
        bool has_negative_distance = false;
		for(int i = 0 ; i < 3 ; i++)
			if(distances[i] > -mEpsilon){
				positive_distances_weight+= mValidRayWeights[i];
                has_positive_distance = true;
            }
			else{
				negative_distances_weight+= mValidRayWeights[i];
                has_negative_distance = true;
            }

        if(has_positive_distance && has_negative_distance)
            mStatistics.NumberOfNoConsensusVotes++;
        else
            mStatistics.NumberOfConsensusVotes++;

		if (positive_distances_weight > negative_distances_weight)
			distance = std::abs(distance);
		else
			distance = -std::abs(distance);

        return distance;
	}

	void CalculateDistanceToSkinProcess::GetRayIntersections(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections)
	{
		//This function passes the ray through the model and gives the hit point to all objects in its way
        //ray is of dimension (3) normalized in (0,1)^3 space
        // direction can be 0,1,2 which are x,y and z respectively

        // first clearing the intersections points vector
        intersections.clear();

        //OctreeType* octree = &mOctree;
        OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess.GetOctreePointer();

        OctreeType::key_type ray_key[3] = {pOctree->CalcKeyNormalized(ray[0]), pOctree->CalcKeyNormalized(ray[1]), pOctree->CalcKeyNormalized(ray[2])};
        OctreeType::key_type cell_key[3];

        // getting the entrance cell from lower extreme
        OctreeType::cell_type* cell = pOctree->pGetCell(ray_key);

        while (cell)
		{
            this->GetCellIntersections(cell, ray, ray_key, direction, intersections);
            // go to the next cell
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key))
			{
                ray_key[direction] = cell_key[direction];
                cell = pOctree->pGetCell(ray_key);
                ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
                //cell get in pGetCell is the right one.
            } else
                cell = NULL;
        }


        // now eliminating the repeated objects
        if (!intersections.empty())
		{
            //sort
            std::sort(intersections.begin(), intersections.end());
            // unique
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (++i_begin != intersections.end())
			{
                // considering the very near points as the same points
                if (fabs(i_begin->first - i_intersection->first) > mEpsilon) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
            intersections.resize((++i_intersection) - intersections.begin());
        }
	}

	int  CalculateDistanceToSkinProcess::GetCellIntersections(OctreeType::cell_type* cell, double* ray,
							 			  					  OctreeType::key_type* ray_key, int direction,
							 		  	  					  std::vector<std::pair<double, Element::GeometryType*> >& intersections)
	{
		//This function passes the ray through the cell and gives the hit point to all objects in its way
		//ray is of dimension (3) normalized in (0,1)^3 space
		// direction can be 0,1,2 which are x,y and z respectively

		typedef OctreeType::cell_type::object_container_type object_container_type;

		object_container_type* objects = (cell->pGetObjects());

		// There are no intersection in empty cells
		if (objects->empty())
			return 0;

		// calculating the two extreme of the ray segment inside the cell
		double ray_point1[3] = {ray[0], ray[1], ray[2]};
		double ray_point2[3] = {ray[0], ray[1], ray[2]};
		double normalized_coordinate;

		// OctreeType* pOctree = (CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess.GetOctree()).get();
		OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess.GetOctreePointer();

		pOctree->CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
		ray_point1[direction] = normalized_coordinate;
		ray_point2[direction] = ray_point1[direction] + pOctree->CalcSizeNormalized(cell);

		pOctree->ScaleBackToOriginalCoordinate(ray_point1);
		pOctree->ScaleBackToOriginalCoordinate(ray_point2);

		for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++)
		{
			double intersection[3]={0.00,0.00,0.00};

			int is_intersected = IntersectionTriangleSegment((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection); // This intersection has to be optimized for axis aligned rays

			if (is_intersected == 1) // There is an intersection but not coplanar
				intersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
		}

		return 0;
	}

	//TODO: This method has been adapted from the previous implementation. It is still pending to update it.
	int CalculateDistanceToSkinProcess::IntersectionTriangleSegment(Element::GeometryType& rGeometry, double* RayPoint1, double* RayPoint2, double* IntersectionPoint)
	{
        array_1d<double,3>    u, v, n;             // triangle vectors
        array_1d<double,3>    dir, w0, w;          // ray vectors
        double     r, a, b;             // params to calc ray-plane intersect


        // get triangle edge vectors and plane normal
        u = rGeometry[1] - rGeometry[0];
        v = rGeometry[2] - rGeometry[0];

        MathUtils<double>::CrossProduct(n, u, v);             // cross product

        if (norm_2(n) == 0)            // triangle is degenerate
            return -1;                 // do not deal with this case

		double triangle_origin_distance = -inner_prod(n, rGeometry[0]);
		Point ray_point_1, ray_point_2;

		for(int i = 0 ; i < 3 ; i++)
        {
            dir[i] = RayPoint2[i] - RayPoint1[i];             // ray direction vector
            w0[i] = RayPoint1[i] - rGeometry[0][i];
			ray_point_1[i] = RayPoint1[i];
			ray_point_2[i] = RayPoint2[i];
		}

		double sign_distance_1 = inner_prod(n, ray_point_1) + triangle_origin_distance;
		double sign_distance_2 = inner_prod(n, ray_point_2) + triangle_origin_distance;

		if (sign_distance_1*sign_distance_2 > mEpsilon) // segment line point on the same side of plane
			return 0;
		a = -inner_prod(n,w0);
        b = inner_prod(n,dir);

        if (fabs(b) < mEpsilon) // ray is parallel to triangle plane
		{
            if (a == 0)                // ray lies in triangle plane
                return 2;
            else return 0;             // ray disjoint from plane
        }

        // get intersect point of ray with triangle plane
        r = a / b;
        if (r < 0.0)                   // ray goes away from triangle
            return 0;                  // => no intersect
        // for a segment, also test if (r > 1.0) => no intersect

        for(int i = 0 ; i < 3 ; i++)
            IntersectionPoint[i]  = RayPoint1[i] + r * dir[i];           // intersect point of ray and plane

        // is I inside T?
        double    uu, uv, vv, wu, wv, D;
        uu = inner_prod(u,u);
        uv = inner_prod(u,v);
        vv = inner_prod(v,v);

        for(int i = 0 ; i < 3 ; i++)
            w[i] = IntersectionPoint[i] - rGeometry[0][i];

        wu = inner_prod(w,u);
        wv = inner_prod(w,v);
        D = uv * uv - uu * vv;

        // get and test parametric coords
        double s, t;
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 - mEpsilon || s > 1.0 + mEpsilon)        // I is outside T
            return 0;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 - mEpsilon || (s + t) > 1.0 + mEpsilon)  // I is outside T
            return 0;

        return 1;                      // I is in T

	}

	void CalculateDistanceToSkinProcess::Execute()
	{
        Internals::Timer timer;
		this->Initialize();
        std::cout << "Initialization is done in " << timer.elapsed().count() << " milliseconds" << std::endl;

        timer.reset();
		this->FindIntersections();
        std::cout << "FindIntersections is done in " << timer.elapsed().count()  << " milliseconds" << std::endl;

        timer.reset();
		this->CalculateDistances(this->GetIntersections());
        std::cout << "CalculateDistances is done in " << timer.elapsed().count()  << " milliseconds" << std::endl;

        		
        ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess).GetModelPart1();

        std::cout << "    NumberOfValidRays        : " << mStatistics.NumberOfValidRays[0] << "," << mStatistics.NumberOfValidRays[1] << "," << mStatistics.NumberOfValidRays[2] << std::endl;
        std::cout << "    NumberOfInvalidRays      : " << mStatistics.NumberOfInvalidRays[0] << "," << mStatistics.NumberOfInvalidRays[1] << "," << mStatistics.NumberOfInvalidRays[2] << std::endl;
        std::cout << "    NumberOfConsensusVotes   : " << mStatistics.NumberOfConsensusVotes << " (" << static_cast<double>(mStatistics.NumberOfConsensusVotes) / ModelPart1.NumberOfNodes() *100.00 << "%)" <<  std::endl;
        std::cout << "    NumberOfNoConsensusVotes : " << mStatistics.NumberOfNoConsensusVotes << " (" << static_cast<double>(mStatistics.NumberOfNoConsensusVotes) / ModelPart1.NumberOfNodes() *100.00 << "%)" << std::endl;
        std::cout << "    NumberOfIrregularElementalDistances : " << mStatistics.NumberOfIrregularElementalDistances << std::endl;

	}

	/// Turn back information as a string.
	std::string CalculateDistanceToSkinProcess::Info() const
	{
		return "CalculateDistanceToSkinProcess";
	}

	/// Print information about this object.
	void CalculateDistanceToSkinProcess::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void CalculateDistanceToSkinProcess::PrintData(std::ostream& rOStream) const
	{
	}



}  // namespace Kratos.
