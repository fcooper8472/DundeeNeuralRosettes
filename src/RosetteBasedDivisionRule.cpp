#include "../../DundeeNeuralRosettes/src/RosetteBasedDivisionRule.hpp"

#include "TypeSixSecretionEnumerations.hpp"
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > RosetteBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    c_vector<double, SPACE_DIM> axis_vector;

    switch (SPACE_DIM)
    {
        case 1:
        {
            EXCEPTION("RosetteBasedDivisionRule is not implemented for SPACE_DIM==1");
        }
        case 2:
        {
        	Node<SPACE_DIM>* p_node = rCellPopulation.GetNodeCorrespondingToCell(pParentCell);

        	const double orientation_angle = p_node->rGetNodeAttributes()[NA_ANGLE];
        	//const double distance = 0.25*p_node->rGetNodeAttributes()[NA_LENGTH] + 0.5*p_node->rGetNodeAttributes()[NA_RADIUS];
        	//const double distance = 0.5*p_node->rGetNodeAttributes()[NA_LENGTH] + p_node->rGetNodeAttributes()[NA_RADIUS];
        	const double distance=0.5*p_node->rGetNodeAttributes()[NA_RADIUS];
            axis_vector(0) = distance*sin(orientation_angle);
            axis_vector(1) = -distance*cos(orientation_angle);
            break;
        }
        case 3:
        {
            EXCEPTION("RosetteBasedDivisionRule is not implemented for SPACE_DIM==3");

        }
        default:
            NEVER_REACHED;
    }

    c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell) - axis_vector;
    c_vector<double, SPACE_DIM> daughter_position = parent_position + 2.0*axis_vector;



    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);
    return positions;
}

// Explicit instantiation
template class RosetteBasedDivisionRule<1,1>;
template class RosetteBasedDivisionRule<1,2>;
template class RosetteBasedDivisionRule<2,2>;
template class RosetteBasedDivisionRule<1,3>;
template class RosetteBasedDivisionRule<2,3>;
template class RosetteBasedDivisionRule<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RosetteBasedDivisionRule)
