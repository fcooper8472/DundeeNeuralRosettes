
#include "../../DundeeNeuralRosettes/src/ApicalForce.hpp"

#include "TypeSixSecretionEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"


#ifdef CHASTE_VTK
#include <vtkLine.h>
#endif // CHASTE_VTK

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <cmath>

#include "Debug.hpp"



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ApicalForce<ELEMENT_DIM, SPACE_DIM>::SetSpringConstant(double springConstant)

{
	mSpringConstant=springConstant;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ApicalForce<ELEMENT_DIM, SPACE_DIM>::GetSpringConstant()

{
	return mSpringConstant;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ApicalForce<ELEMENT_DIM, SPACE_DIM>::SetSpringLength(double springLength)

{
	mSpringLength=springLength;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ApicalForce<ELEMENT_DIM, SPACE_DIM>::GetSpringLength()

{
	return mSpringLength;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ApicalForce<ELEMENT_DIM, SPACE_DIM>::ApicalForce()
        : AbstractForce<ELEMENT_DIM, SPACE_DIM>(),
          mSpringConstant(10.0), mSpringLength(1.0)
{
    assert(ELEMENT_DIM == 2u);
    assert(SPACE_DIM == 2u);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ApicalForce<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    auto p_cell_population = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);

    if (p_cell_population == nullptr)
    {
        EXCEPTION("Capsule force only works with AbstractCentreBasedCellPopulation");
    }
    double 	R_rosette= 20.0;

    // Set all applied angles back to zero
    for (auto iter = p_cell_population->rGetMesh().GetNodeIteratorBegin();
         iter != p_cell_population->rGetMesh().GetNodeIteratorEnd();
         ++iter)
    {


    	unsigned node_index = iter->GetIndex();

        // resting spring length is function of cell cycle state

    	auto p_cell_cycle_model=static_cast<AbstractPhaseBasedCellCycleModel*>((p_cell_population->GetCellUsingLocationIndex(node_index))->GetCellCycleModel());


        if (p_cell_cycle_model== nullptr)
        {
        	mSpringLength=4.0;
        }
        else
        {

        switch (p_cell_cycle_model->GetCurrentCellCyclePhase())
		{

			case G_ZERO_PHASE:
				mSpringLength=15.0;
				break;
			case G_ONE_PHASE:
				mSpringLength=15.0;
				//PRINT_VARIABLE(p_cell_cycle_model->GetCurrentCellCyclePhase());

				break;

			case S_PHASE:
				mSpringLength=4.5;
				break;

			case G_TWO_PHASE:
				mSpringLength=4.5;
				break;


			case M_PHASE:
				mSpringLength=15.0;
				break;


			default:


				NEVER_REACHED;
		}
        }







            c_vector<double, SPACE_DIM> nodelocation=iter->rGetLocation();

            const double angle_a = iter->rGetNodeAttributes()[NA_ANGLE];
            const double radius = iter->rGetNodeAttributes()[NA_RADIUS];
            const double length_a = iter->rGetNodeAttributes()[NA_LENGTH]+2.0*radius;


            c_vector<double, SPACE_DIM> apicalendlocation;

            apicalendlocation[0]=  nodelocation[0] + 0.5 * length_a * cos(angle_a);
            apicalendlocation[1]=  nodelocation[1] + 0.5 * length_a * sin(angle_a);

            c_vector<double, SPACE_DIM> apical_foot_of_cell;
            apical_foot_of_cell[0]=0;
            apical_foot_of_cell[1]=0;





            double force_magnitude = mSpringConstant*(norm_2(apicalendlocation-apical_foot_of_cell)-mSpringLength);

            c_vector<double, SPACE_DIM> force_direction_a_to_b =apicalendlocation/norm_2(apicalendlocation) ;

            c_vector<double, SPACE_DIM> force_a_b = -force_direction_a_to_b * force_magnitude;


            // Calculate the 2D cross product of two vectors
            auto cross_product = [](c_vector<double, SPACE_DIM> a, c_vector<double, SPACE_DIM> b) -> double
            {
                return a[0] * b[1] - b[0] * a[1];
            };

            iter->rGetNodeAttributes()[NA_APPLIED_ANGLE] += cross_product(-apicalendlocation+nodelocation, force_a_b);

            iter->AddAppliedForceContribution(force_a_b);


            c_vector<double, SPACE_DIM> basalendlocation;

            basalendlocation[0]=  nodelocation[0] - 0.5 * length_a * cos(angle_a);
            basalendlocation[1]=  nodelocation[1] - 0.5 * length_a * sin(angle_a);

            double cell_angle_wrt_origin=atan2(nodelocation[1],nodelocation[0] );

                       c_vector<double, SPACE_DIM> basal_foot_of_cell;
                       basal_foot_of_cell[0]=R_rosette*cos(cell_angle_wrt_origin);
                       basal_foot_of_cell[1]=R_rosette*sin(cell_angle_wrt_origin);

                       double basal_SpringLength= R_rosette-mSpringLength;





                       double basal_force_magnitude = mSpringConstant*(norm_2(basalendlocation-basal_foot_of_cell)-basal_SpringLength);

                       c_vector<double, SPACE_DIM> basal_force_direction_a_to_b =(basalendlocation-basal_foot_of_cell)/norm_2(basalendlocation-basal_foot_of_cell) ;

                       c_vector<double, SPACE_DIM> basal_force_a_b = -basal_force_magnitude * basal_force_direction_a_to_b;




                       iter->rGetNodeAttributes()[NA_APPLIED_ANGLE] += cross_product(-basalendlocation+nodelocation, basal_force_a_b);

                       iter->AddAppliedForceContribution(basal_force_a_b);






    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ApicalForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ApicalForce<1,1>;
template class ApicalForce<1,2>;
template class ApicalForce<2,2>;
template class ApicalForce<1,3>;
template class ApicalForce<2,3>;
template class ApicalForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ApicalForce)
