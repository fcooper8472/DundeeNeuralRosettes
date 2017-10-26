
#ifndef TESTROSETTESIMULATION_HPP_
#define TESTROSETTESIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <cycle/UniformCellCycleModel.hpp>
#include "../../DundeeNeuralRosettes/src/ApicalForce.hpp"
#include "../../DundeeNeuralRosettes/src/NodeBasedCellPopulationWithCapsulesRosettes.hpp"
#include "../../DundeeNeuralRosettes/src/RosetteBasedDivisionRule.hpp"

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "NoCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Cell.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "UniformCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellIdWriter.hpp"

// Header files included in this project
#include "TypeSixSecretionEnumerations.hpp"
#include "ForwardEulerNumericalMethodForCapsulesRosettes.hpp"
#include "CapsuleForce.hpp"
#include "DiffusionForce.hpp"



#include "CapsuleOrientationWriter.hpp"
#include "CapsuleScalingWriter.hpp"
#include "SquareBoundaryCondition.hpp"
#include "CapsuleBasedDivisionRule.hpp"
#include "TypeSixMachineModifier.hpp"
#include "NodeBasedCellPopulationWithCapsules.hpp"
#include "TypeSixMachineProperty.hpp"
#include "TypeSixMachineCellKiller.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GammaG1CellCycleModel.hpp"



#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestRosetteSimulation : public AbstractCellBasedTestSuite
{
public:






    void NoTestRosettes() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        const unsigned num_nodes = 8u;
        auto p_rand_gen = RandomNumberGenerator::Instance();

        // Create some capsules
        std::vector<Node<2>*> nodes;

        double R_rosettes=5.0;
        for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
        {
            c_vector<double, 2> safe_location;
            double node_angle=  2.0*M_PI*double(node_idx)/double(num_nodes)-M_PI;
            double cell_centre_rad=R_rosettes;
            safe_location = Create_c_vector(cell_centre_rad*cos(node_angle), cell_centre_rad*sin(node_angle));
            nodes.push_back(new Node<2>(node_idx, safe_location));
        }

        /*
         * We then convert this list of nodes to a `NodesOnlyMesh`,
         * which doesn't do very much apart from keep track of the nodes.
         */
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 10.0);



        for (unsigned node_idx = 0; node_idx < mesh.GetNumNodes(); ++node_idx)
        {
            mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
            mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
            double node_angle=  2.0*M_PI*double(node_idx)/double(num_nodes)-M_PI;

            double rand_angle = 0.3*(RandomNumberGenerator::Instance()->ranf()- 0.5);

            mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_ANGLE] = node_angle+rand_angle;
            mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = p_rand_gen->NormalRandomDeviate(2.0, 0.5);
            mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
        }

        //Create cells
        std::vector<CellPtr> cells;
        auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        // Create cell population
        NodeBasedCellPopulation<2> population(mesh, cells);

        population.AddCellWriter<CapsuleOrientationWriter>();
        population.AddCellWriter<CapsuleScalingWriter>();

        // Create simulation
        OffLatticeSimulation<2> simulator(population);
        simulator.SetOutputDirectory("TestRosetteInitialisation");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(1u);

        auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsulesRosettes<2,2>>();
        simulator.SetNumericalMethod(p_numerical_method);

        /*
         * We now create a force law and pass it to the simulation
         * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
         */
        auto p_calsule_force = boost::make_shared<CapsuleForce<2>>();
        simulator.AddForce(p_calsule_force);

        auto p_apical_force = boost::make_shared<ApicalForce<2>>();
        simulator.AddForce(p_apical_force);

      //simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /* We then set an end time and run the simulation */
        simulator.SetEndTime(5.0);
        simulator.Solve();
    }

    void TestRosettesWithDivision() throw (Exception)
        {
            EXIT_IF_PARALLEL;

            const unsigned num_nodes = 8u;
            //auto p_rand_gen = RandomNumberGenerator::Instance();

            // Create some capsules
            std::vector<Node<2>*> nodes;

            double R_rosettes=5.0;
            for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
            {
                c_vector<double, 2> safe_location;
                double node_angle=  2.0*M_PI*double(node_idx)/double(num_nodes)-M_PI;
                double cell_centre_rad=R_rosettes;
                safe_location = Create_c_vector(cell_centre_rad*cos(node_angle), cell_centre_rad*sin(node_angle));
                nodes.push_back(new Node<2>(node_idx, safe_location));
            }

            /*
             * We then convert this list of nodes to a `NodesOnlyMesh`,
             * which doesn't do very much apart from keep track of the nodes.
             */
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 10.0);



            for (unsigned node_idx = 0; node_idx < mesh.GetNumNodes(); ++node_idx)
            {
                mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
                mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
                double node_angle=  2.0*M_PI*double(node_idx)/double(num_nodes)-M_PI;

                double rand_angle = 0.3*(RandomNumberGenerator::Instance()->ranf()- 0.5);

                mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_ANGLE] = node_angle+rand_angle;
                mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
                mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 1.0;
            }



            std::vector<CellPtr> cells;
           		  MAKE_PTR(WildTypeCellMutationState, p_state);
           		  MAKE_PTR(TransitCellProliferativeType, p_type);
           		  for (unsigned i=0; i<mesh.GetNumNodes(); i++)
           		  {
           			 //FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
           			 //UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();
           			GammaG1CellCycleModel* p_model = new GammaG1CellCycleModel();

           			 p_model->SetShape(4.517);
           			 p_model->SetScale(1.986);

           			CellPtr p_cell(new Cell(p_state, p_model));
           			  p_cell->SetCellProliferativeType(p_type);


           			  p_model->SetG1Duration();

           			  p_model->SetSDuration(2.0);
           			  p_model->SetG2Duration(2.0);
           			  p_model->SetMDuration(0.5);

           			  double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_model->GetG1Duration()+p_model->GetSG2MDuration() );
           			  PRINT_VARIABLE(birth_time);
           			  p_cell->SetBirthTime(birth_time);

					  double length = mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH];
					  double initial_radius=1.0;

					  double division_radius = -length/M_PI+sqrt(pow(length,2)/pow(M_PI,2)+2.0*pow(initial_radius,2)+4.0*initial_radius*length/M_PI); //% + 2*pNodeA->rGetNodeAttributes()[NA_RADIUS];

           			  mesh.GetNode(i)->rGetNodeAttributes()[NA_RADIUS] = initial_radius +(division_radius-initial_radius)*p_cell->GetAge()/p_model->GetAverageTransitCellCycleTime(); ;





           			  cells.push_back(p_cell);
           		  }






            // Create cell population
            NodeBasedCellPopulationWithCapsulesRosettes<2> population(mesh, cells);

            population.AddCellWriter<CapsuleOrientationWriter>();
            population.AddCellWriter<CapsuleScalingWriter>();

            // Create simulation
            OffLatticeSimulation<2> simulator(population);
            simulator.SetOutputDirectory("TestRosetteWithDivision");
            simulator.SetDt(1.0/520.0);
            simulator.SetSamplingTimestepMultiple(20u);

            auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsulesRosettes<2,2>>();
            p_numerical_method->SetAxialCapsuleGrowth(false);
            simulator.SetNumericalMethod(p_numerical_method);

            boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new RosetteBasedDivisionRule<2,2>());
            population.SetCentreBasedDivisionRule(p_division_rule);


            /*
             * We now create a force law and pass it to the simulation
             * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
             */
            auto p_calsule_force = boost::make_shared<CapsuleForce<2>>();
            simulator.AddForce(p_calsule_force);

            auto p_apical_force = boost::make_shared<ApicalForce<2>>();
            p_apical_force->SetSpringConstant(20.01);
            simulator.AddForce(p_apical_force);

           auto p_diffusion_force = boost::make_shared<DiffusionForce<2>>();
            p_diffusion_force->SetAbsoluteTemperature(290.30);
           simulator.AddForce(p_diffusion_force);

          //simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

            /* We then set an end time and run the simulation */
            simulator.SetEndTime(38.0);

            // Check the correct solution was obtained
            //        for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
              //           cell_iter != simulator.rGetCellPopulation().End();
                //         ++cell_iter)
                  //  {
                        // Get cell model
                    //    AbstractCellCycleModel* p_abstract_model = cell_iter->GetCellCycleModel();
                      //  FixedG1GenerationalCellCycleModel* p_model = static_cast<FixedG1GenerationalCellCycleModel*> (p_abstract_model);
                        //PRINT_VARIABLE(p_model->GetAverageTransitCellCycleTime());
                        //PRINT_VARIABLE(p_model->GetCurrentCellCyclePhase());
                    //}

            simulator.Solve();

            //PRINT_VARIABLE(simulator.rGetPopulation().GetNumRealCells());
        }


};

#endif /*TESTROSETTESIMULATION_HPP_*/
