// Run an EMI model simulation of a cerebellar purkinje neuron

#include "mfem.hpp"
#include <omp.h>
#include <fstream>
#include <iostream>
#include "Masoli2015.h"
#include <chrono>

using namespace std;
using namespace mfem;

int load_number(std::string load_string);
void load_number_cells(std::string load_string, int *some_number, int num_cells);
void load_array(std::string load_string, Array<int> &some_list, int list_length);
void load_array_cells(std::string load_string, Array<int> *some_list, int *list_length, int num_cells);
void load_integration_rule(std::string load_string, IntegrationRule &some_list, int list_length);
void load_integration_rule_cells(std::string load_string, IntegrationRule *some_list, int *list_length, int num_cells);
void write_mean_values(std::ofstream &wf, Vector sol_vec, int num_cells, int* num_mem_nodes_cell, Array<int> NU_nodes, Array<int> *membrane_list_cell, FiniteElementSpace **fespace_cell_serial, Array<int> mem_attr);


int main(int argc, char *argv[])
{
    // Initialize MPI and HYPRE.
    Mpi::Init();
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
#if defined(USE_HYPRE) || defined(USE_SUPERLU)
    Hypre::Init();
#endif
    
    // Parse command-line options.
    double Tstop = 500;          // Total simulation time (in ms)
    double dt = 0.01;            // Global time step (in ms)
    double dt_ode = 1.0e-4;      // Membrane model time step (in ms)
    bool save_paraview = false;  // Save solution to ParaView
    bool U_dist = false;         // Uniform ion channel distribution
    
    OptionsParser args(argc, argv);
    args.AddOption(&save_paraview, "-s", "--save_paraview", "-no-s", "--no_save_paraview", "Save paraview solution");
    args.AddOption(&U_dist, "-u", "--U_dist", "-no-u", "--no_U_dist", "Uniform channel distribution");
    args.AddOption(&Tstop, "-t", "--Tstop", "Simulation time.");
    args.AddOption(&dt, "-dt", "--dt", "Time step.");
    
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    if (myid == 0) {
        args.PrintOptions(cout);
    }
    
    // Set up simulation title
    std::string simulation_name = "one_cell";
    if (U_dist) {
        simulation_name += "_uniform";
    }
    
    // Model parameter values
    const double Cm = 1.0;
    const double sigma_e = 3.0;
    const double sigma_i = 8.2;

    // Specify time stepping
    int save_limit = 1;
    int print_limit = 10;
    int Nt = (int)round(Tstop/dt);
    int M_it = 1;  // Operator splitting inner iterations
    int N_it = 1;  // Operator splitting outer iterations
    if (dt_ode > dt) {
        dt_ode = dt;
    }
    int K = (int)round(dt/dt_ode);
    
    // Specify mesh
    int num_cells = 1;
    int dim = 3;
    const char *intracellular_mesh_file = "meshes/one_cell/intracellular.msh";
    const char *extracellular_mesh_file = "meshes/one_cell/extracellular.msh";
    std::string cell_base_file = "meshes/one_cell/cell";
    std::string load_mesh_string = "meshes/mesh_properties/one_cell_7branch";
    std::string current_load_string;
    int current_idx;
    
    if (U_dist) {
        load_mesh_string += "_no_mye";
    }
    
    // Load meshes
    Mesh *mesh_e = new Mesh(extracellular_mesh_file, 1, 1);
    Mesh *mesh_i = new Mesh(intracellular_mesh_file, 1, 1);
    std::string cell_string;
    Mesh *mesh_cell[num_cells];
    int* cell_partitioning[num_cells];
    for (int i = 0; i < num_cells; i++)
    {
        cell_string.assign(cell_base_file);
        cell_string += to_string(i+1);
        cell_string += ".msh";
        mesh_cell[i] = new Mesh(cell_string.c_str(), 1, 1);
        cell_partitioning[i] = mesh_cell[i]->GeneratePartitioning(num_procs);
    }
    ParMesh *pmesh_cell[num_cells];
    for (int i = 0; i < num_cells; i++)
    {
        pmesh_cell[i] = new ParMesh(MPI_COMM_WORLD, *mesh_cell[i], cell_partitioning[i]);
    }
    ParMesh *pmesh_e = new ParMesh(MPI_COMM_WORLD, *mesh_e);
    
    
    
    // Enable hardware devices
    const char *device_config = "cpu";
    Device device(device_config);
    
    // Define a finite element space on the mesh
    int order = 1;
    FiniteElementCollection *fec_e;
    FiniteElementCollection *fec_i;
    FiniteElementCollection *fec_cell[num_cells];
    fec_e = new H1_FECollection(order, dim);
    fec_i = new H1_FECollection(order, dim);
    ParFiniteElementSpace *fespace_e = new ParFiniteElementSpace(pmesh_e, fec_e);
    ParFiniteElementSpace *fespace_cell[num_cells];
    FiniteElementSpace *fespace_cell_serial[num_cells];
    MPI_Comm comm = pmesh_e->GetComm();
    FiniteElementSpace *fespace_i = new FiniteElementSpace(mesh_i, fec_i);
    if (myid == 0) {
        cout << "Number of finite element unknowns (extracellular): "
        << fespace_e->GetTrueVSize() << endl;
        cout << "Number of finite element unknowns (intracellular): "
        << fespace_i->GetTrueVSize() << endl;
    }
    
    
    for (int i = 0; i < num_cells; i++)
    {
        fec_cell[i] = new H1_FECollection(order, dim);
        fespace_cell[i] = new ParFiniteElementSpace(pmesh_cell[i], fec_cell[i]);
        fespace_cell_serial[i] = new FiniteElementSpace(mesh_cell[i], fec_cell[i]);
    }
    
    

   // Set up mapping between membrane points in different meshes
   current_load_string = load_mesh_string;
   current_load_string += "_num_mem_nodes_i";
   int num_mem_nodes_i = load_number(current_load_string);
   current_load_string = load_mesh_string;
   current_load_string += "_num_mem_nodes_e";
   int num_mem_nodes_e = load_number(current_load_string);
   Array<int> membrane_list_e;
   Array<int> membrane_elements_e;
   IntegrationRule membrane_rule_e(num_mem_nodes_e);
   current_load_string = load_mesh_string;
   current_load_string += "_membrane_list_e";
   load_array(current_load_string, membrane_list_e, num_mem_nodes_e);
   current_load_string = load_mesh_string;
   current_load_string += "_membrane_elements_e";
   load_array(current_load_string, membrane_elements_e, num_mem_nodes_e);
   current_load_string = load_mesh_string;
   current_load_string += "_membrane_rule_e";
   load_integration_rule(current_load_string, membrane_rule_e, num_mem_nodes_e);
    
    if (myid == 0) {
        cout << "Number of finite element unknowns (membrane): "
        << num_mem_nodes_i << endl;
    }
    
    // Set up mapping for the membrane of the cells
    int num_mem_nodes_cell[num_cells];
    Array<int> membrane_list_cell[num_cells];
    Array<int> membrane_idx_for_cell[num_cells];
    Array<int> membrane_elements_cell[num_cells];
    IntegrationRule membrane_rule_cell[num_cells];
    
    current_load_string = load_mesh_string;
    current_load_string += "_num_mem_nodes_total";
    int num_mem_nodes_total = load_number(current_load_string);
    current_load_string = load_mesh_string;
    current_load_string += "_num_mem_nodes_cell";
    load_number_cells(current_load_string, num_mem_nodes_cell, num_cells);
    current_load_string = load_mesh_string;
    current_load_string += "_membrane_list_cell";
    load_array_cells(current_load_string, membrane_list_cell, num_mem_nodes_cell, num_cells);
    current_load_string = load_mesh_string;
    current_load_string += "_membrane_idx_for_cell";
    load_array_cells(current_load_string, membrane_idx_for_cell, num_mem_nodes_cell, num_cells);
    
    current_load_string = load_mesh_string;
    current_load_string += "_membrane_elements_cell";
    load_array_cells(current_load_string, membrane_elements_cell, num_mem_nodes_cell, num_cells);
    current_load_string = load_mesh_string;
    current_load_string += "_membrane_rule_cell";
    load_integration_rule_cells(current_load_string, membrane_rule_cell, num_mem_nodes_cell, num_cells);
    
    
    // Set up mapping from a cell boundary to the intracellular space
    Array<int> cell_for_i_point;
    Array<int> element_for_i_point;
    IntegrationRule rule_for_i_point;
    current_load_string = load_mesh_string;
    current_load_string += "_num_nodes_i";
    int num_nodes_i = load_number(current_load_string);
    current_load_string = load_mesh_string;
    current_load_string += "_cell_for_i_point";
    load_array(current_load_string, cell_for_i_point, num_nodes_i);
    current_load_string = load_mesh_string;
    current_load_string += "_element_for_i_point";
    load_array(current_load_string, element_for_i_point, num_nodes_i);
    current_load_string = load_mesh_string;
    current_load_string += "_rule_for_i_point";
    load_integration_rule(current_load_string, rule_for_i_point, num_nodes_i);
    
    // Load NU nodes
    Array<int> NU_nodes;
    current_load_string = load_mesh_string;
    current_load_string += "_NU_nodes";
    load_array(current_load_string, NU_nodes, num_mem_nodes_total);

   // Determine the list of true (i.e. conforming) essential boundary dofs
   Array<int> ess_tdof_list_e;
   Array<int> ess_tdof_list_cell[num_cells];
   Array<int> ess_bdr_cell[num_cells];
   Array<int> ess_bdr_e(pmesh_e->bdr_attributes.Max());

   // Extracellular Dirichlet boundary condition on nodes marked by attribute 1
   ess_bdr_e = 0;
   ess_bdr_e[0] = 1;
   fespace_e->GetEssentialTrueDofs(ess_bdr_e, ess_tdof_list_e);

   // No intracellular Dirichlet boundary conditions
   for (int i = 0; i < num_cells; i++)
   {
       ess_bdr_cell[i].SetSize(pmesh_cell[i]->bdr_attributes.Max());
       ess_bdr_cell[i] = 0;
       fespace_cell[i]->GetEssentialTrueDofs(ess_bdr_cell[i], ess_tdof_list_cell[i]);
   }

   // Define the solution vectors
   GridFunction x_i(fespace_i);
   ParGridFunction x_cell[num_cells];
   ParGridFunction x_e(fespace_e);
   GridFunction v_in_i(fespace_i);
    
   //  Initial conditions
   x_e = 0.0;
   x_i = 0.0;
   v_in_i = 0.0;

   for (int i = 0; i < num_cells; i++)
   {
       x_cell[i].SetSpace(fespace_cell[i]);
       x_cell[i] = 0.0;
   }
    
    GridFunction x_cell_serial[num_cells];
    for (int i=0; i < num_cells; i++) {
        x_cell_serial[i].SetSpace(fespace_cell_serial[i]);
    }
    
    // Membrane attribute arrays
    Array<int> mem_attr;
    mem_attr.SetSize(mesh_i->bdr_attributes.Max());
    mem_attr = 1;
    mem_attr[0] = 0; // Myelin marked 3
    mem_attr[1] = 0; // Myelin marked 3
    if (U_dist == false) {
        mem_attr[2] = 0; // Myelin marked 3
    }



    // Coefficients for bilinear forms
    ConstantCoefficient Cmdt_coeff(Cm/dt);
    ConstantCoefficient sigma_i_coeff(sigma_i);
    ConstantCoefficient sigma_e_coeff(sigma_e);

    // Set up functions and vectors for rhs expressions
    ParGridFunction b_e_func(fespace_e);
    GridFunctionCoefficient b_e_coeff(&b_e_func);
    Vector b_e_vec(fespace_e->GetTrueVSize());
    Vector b_e_mem_vec(num_mem_nodes_e);
    LinearForm *b_e = new LinearForm(fespace_e);
    b_e->AddBoundaryIntegrator(new BoundaryLFIntegrator(b_e_coeff), mem_attr);

    Vector b_i_vec(fespace_i->GetTrueVSize());
    Vector b_i_mem_vec(num_mem_nodes_i);

    Vector ui_prev_vec_e(num_mem_nodes_e);
    Vector ui_prev_vec_cell(num_mem_nodes_total);
    Vector v_prev_vec(num_mem_nodes_total);
    Vector v_prev_vec_e(num_mem_nodes_e);
    Vector ue_prev_vec_e(num_mem_nodes_e);
    Vector ue_prev_vec_cell(num_mem_nodes_total);
    Vector ue_new_vec_cell(num_mem_nodes_total);

    LinearForm *b_cell[num_cells];
    ParGridFunction *b_cell_func[num_cells];
    GridFunction b_cell_func_serial[num_cells];
    GridFunction b_w_cell_func[num_cells];
    GridFunction v_in_cell[num_cells];
    Vector b_cell_vec[num_cells];
    Vector b_cell_bdr_vec[num_cells];
    Vector x_cell_vec[num_cells];
    

    for (int i = 0; i < num_cells; i++)
    {
        b_cell_func_serial[i].SetSpace(fespace_cell_serial[i]);
        v_in_cell[i].SetSpace(fespace_cell_serial[i]);
        b_cell_vec[i].SetSize(fespace_cell[i]->GetTrueVSize());
        x_cell[i].GetTrueDofs(x_cell_vec[i]);
    }
    

    // Update initial conditions
    ue_prev_vec_e = 0.0;
    ue_prev_vec_cell = 0.0;
    ue_new_vec_cell = 0.0;

    
    // Membrane ODE model (Masoli2015)
    const int num_states = 46;
    const int num_param = 28;
    int V_idx = 45;
    int Na_idx = 9;
    int Kv11_idx = 11;
    int Kv15_idx = 12;
    int Kv33_idx = 13;
    int Kv34_idx = 14;
    int Kv43_idx = 15;
    int Kir2x_idx = 16;
    int Kca11_idx = 17;
    int Kca22_idx = 18;
    int Kca31_idx = 19;
    int Cav21_idx = 20;
    int Cav31_idx = 21;
    int Cav32_idx = 22;
    int Cav33_idx = 23;
    int HCN1_idx = 25;
    int leak_idx = 27;
    int Tdiff_idx = 4;
    
    // Initialize ODE model parameters
    double parameters[num_param];
    double (*parameters_nodes)[num_param] = new double[num_mem_nodes_total][num_param];
    init_parameters_values(parameters);
    double init_states[num_states];
    init_state_values(init_states);


    // Set up state variables
    double (*states)[num_states] = new double[num_mem_nodes_total][num_states];
    current_idx = 0;
    for (int i = 0; i < num_cells; i++)
    {
        for (int k=0; k<num_mem_nodes_cell[i]; k++) {
            for (int j = 0; j < num_states; j++)
            {
                states[current_idx][j] = init_states[j];
            }
            
            for (int j = 0; j < num_param; j++)
            {
                parameters_nodes[current_idx][j] = parameters[j];
            }
            
            // Specific parameters
            if (U_dist == false) {
                // Dendrite r > 4 um
                if (NU_nodes[current_idx] == 4) {
                    parameters_nodes[current_idx][Na_idx] = 0.016;
                    parameters_nodes[current_idx][Kv11_idx] = 0.0012;
                    parameters_nodes[current_idx][Kv15_idx] = 1.3e-4;
                    parameters_nodes[current_idx][Kv33_idx] = 0.01;
                    parameters_nodes[current_idx][Kv34_idx] = 0;
                    parameters_nodes[current_idx][Kv43_idx] = 0.001;
                    parameters_nodes[current_idx][Kir2x_idx] = 0.00001;
                    parameters_nodes[current_idx][Kca11_idx] = 3.5e-2;
                    parameters_nodes[current_idx][Kca22_idx] = 1.0e-3;
                    parameters_nodes[current_idx][Kca31_idx] = 0.002;
                    parameters_nodes[current_idx][Cav21_idx] = 1.0e-3;
                    parameters_nodes[current_idx][Cav31_idx] = 5.0e-6;
                    parameters_nodes[current_idx][Cav32_idx] = 0.0012;
                    parameters_nodes[current_idx][Cav33_idx] = 0.0001;
                    parameters_nodes[current_idx][HCN1_idx] = 0.000004;
                    parameters_nodes[current_idx][leak_idx] = 0.0;
                    
                    // Dendrite r > 1.75 um
                } else if (NU_nodes[current_idx] == 5) {
                    parameters_nodes[current_idx][Na_idx] = 0.0;
                    parameters_nodes[current_idx][Kv11_idx] = 0.0012;
                    parameters_nodes[current_idx][Kv15_idx] = 1.3e-4;
                    parameters_nodes[current_idx][Kv33_idx] = 0.01;
                    parameters_nodes[current_idx][Kv34_idx] = 0;
                    parameters_nodes[current_idx][Kv43_idx] = 0.001;
                    parameters_nodes[current_idx][Kir2x_idx] = 0.00001;
                    parameters_nodes[current_idx][Kca11_idx] = 3.5e-2;
                    parameters_nodes[current_idx][Kca22_idx] = 1.0e-3;
                    parameters_nodes[current_idx][Kca31_idx] = 0.002;
                    parameters_nodes[current_idx][Cav21_idx] = 1.0e-3;
                    parameters_nodes[current_idx][Cav31_idx] = 5.0e-6;
                    parameters_nodes[current_idx][Cav32_idx] = 0.0012;
                    parameters_nodes[current_idx][Cav33_idx] = 0.0001;
                    parameters_nodes[current_idx][HCN1_idx] = 0.000004;
                    parameters_nodes[current_idx][leak_idx] = 0.0;
                    
                    // Dendrite r < 1.75 um
                } else if (NU_nodes[current_idx] == 6) {
                    parameters_nodes[current_idx][Na_idx] = 0.0;
                    parameters_nodes[current_idx][Kv11_idx] = 0.0012*30;
                    parameters_nodes[current_idx][Kv15_idx] = 1.3e-4*30;
                    parameters_nodes[current_idx][Kv33_idx] = 0.01*30;
                    parameters_nodes[current_idx][Kv34_idx] = 0;
                    parameters_nodes[current_idx][Kv43_idx] = 0.001*30;
                    parameters_nodes[current_idx][Kir2x_idx] = 0.0;
                    parameters_nodes[current_idx][Kca11_idx] = 3.5e-2*30;
                    parameters_nodes[current_idx][Kca22_idx] = 0.0;
                    parameters_nodes[current_idx][Kca31_idx] = 0.0;
                    parameters_nodes[current_idx][Cav21_idx] = 1.0e-3;
                    parameters_nodes[current_idx][Cav31_idx] = 0.0;
                    parameters_nodes[current_idx][Cav32_idx] = 0.0;
                    parameters_nodes[current_idx][Cav33_idx] = 0.0001;
                    parameters_nodes[current_idx][HCN1_idx] = 0.000004;
                    parameters_nodes[current_idx][leak_idx] = 0.0;
                    
                    // AIS
                } else if (NU_nodes[current_idx] == 8) {
                    parameters_nodes[current_idx][Na_idx] = 1.5;
                    parameters_nodes[current_idx][Kv11_idx] = 0.0;
                    parameters_nodes[current_idx][Kv15_idx] = 0.0;
                    parameters_nodes[current_idx][Kv33_idx] = 0.0;
                    parameters_nodes[current_idx][Kv34_idx] = 0.01;
                    parameters_nodes[current_idx][Kv43_idx] = 0.0;
                    parameters_nodes[current_idx][Kir2x_idx] = 0.0;
                    parameters_nodes[current_idx][Kca11_idx] = 0.0;
                    parameters_nodes[current_idx][Kca22_idx] = 0.0;
                    parameters_nodes[current_idx][Kca31_idx] = 0.0;
                    parameters_nodes[current_idx][Cav21_idx] = 2.2e-4;
                    parameters_nodes[current_idx][Cav31_idx] = 1.0e-5;
                    parameters_nodes[current_idx][Cav32_idx] = 0.0;
                    parameters_nodes[current_idx][Cav33_idx] = 0.0;
                    parameters_nodes[current_idx][HCN1_idx] = 0.0;
                    parameters_nodes[current_idx][leak_idx] = 0.0;
                
                // Para AIS
                } else if (NU_nodes[current_idx] == 9) {
                    parameters_nodes[current_idx][Na_idx] = 0.0;
                    parameters_nodes[current_idx][Kv11_idx] = 0.01;
                    parameters_nodes[current_idx][Kv15_idx] = 0.0;
                    parameters_nodes[current_idx][Kv33_idx] = 0.0;
                    parameters_nodes[current_idx][Kv34_idx] = 0.0;
                    parameters_nodes[current_idx][Kv43_idx] = 0.0;
                    parameters_nodes[current_idx][Kir2x_idx] = 0.0;
                    parameters_nodes[current_idx][Kca11_idx] = 0.0;
                    parameters_nodes[current_idx][Kca22_idx] = 0.0;
                    parameters_nodes[current_idx][Kca31_idx] = 0.0;
                    parameters_nodes[current_idx][Cav21_idx] = 0.0;
                    parameters_nodes[current_idx][Cav31_idx] = 0.0;
                    parameters_nodes[current_idx][Cav32_idx] = 0.0;
                    parameters_nodes[current_idx][Cav33_idx] = 0.0;
                    parameters_nodes[current_idx][HCN1_idx] = 0.0;
                    parameters_nodes[current_idx][leak_idx] = 0.0;
                
                // Ranvier nodes
                } else if (NU_nodes[current_idx] == 10) {
                    parameters_nodes[current_idx][Na_idx] = 0.03;
                    parameters_nodes[current_idx][Kv11_idx] = 0.0;
                    parameters_nodes[current_idx][Kv15_idx] = 0.0;
                    parameters_nodes[current_idx][Kv33_idx] = 0.0;
                    parameters_nodes[current_idx][Kv34_idx] = 0.01;
                    parameters_nodes[current_idx][Kv43_idx] = 0.0;
                    parameters_nodes[current_idx][Kir2x_idx] = 0.0;
                    parameters_nodes[current_idx][Kca11_idx] = 0.0;
                    parameters_nodes[current_idx][Kca22_idx] = 0.0;
                    parameters_nodes[current_idx][Kca31_idx] = 0.0;
                    parameters_nodes[current_idx][Cav21_idx] = 2.2e-4;
                    parameters_nodes[current_idx][Cav31_idx] = 1.0e-5;
                    parameters_nodes[current_idx][Cav32_idx] = 0.0;
                    parameters_nodes[current_idx][Cav33_idx] = 0.0;
                    parameters_nodes[current_idx][HCN1_idx] = 0.0;
                    parameters_nodes[current_idx][leak_idx] = 0.0;
                
                // Collateral
                } else if (NU_nodes[current_idx] == 11) {
                    parameters_nodes[current_idx][Na_idx] = 0.03;
                    parameters_nodes[current_idx][Kv11_idx] = 0.0;
                    parameters_nodes[current_idx][Kv15_idx] = 0.0;
                    parameters_nodes[current_idx][Kv33_idx] = 0.0;
                    parameters_nodes[current_idx][Kv34_idx] = 0.02;
                    parameters_nodes[current_idx][Kv43_idx] = 0.0;
                    parameters_nodes[current_idx][Kir2x_idx] = 0.0;
                    parameters_nodes[current_idx][Kca11_idx] = 0.0;
                    parameters_nodes[current_idx][Kca22_idx] = 0.0;
                    parameters_nodes[current_idx][Kca31_idx] = 0.0;
                    parameters_nodes[current_idx][Cav21_idx] = 2.2e-4;
                    parameters_nodes[current_idx][Cav31_idx] = 1.0e-5;
                    parameters_nodes[current_idx][Cav32_idx] = 0.0;
                    parameters_nodes[current_idx][Cav33_idx] = 0.0;
                    parameters_nodes[current_idx][HCN1_idx] = 0.0;
                    parameters_nodes[current_idx][leak_idx] = 0.0;
                }
                
            }
            
            current_idx++;
        }
    }
    

    // Update initial conditions
    current_idx = 0;
    for (int i=0; i < num_cells; i++) {
        x_cell_serial[i] = 0.0;
        for (int k=0; k<num_mem_nodes_cell[i]; k++) {
            x_cell_serial[i][membrane_list_cell[i][k]] = states[current_idx][V_idx];
            current_idx++;
        }
    }
    
    for (int i = 0; i < fespace_i->GetNDofs(); i++)
    {
        x_i[i] = x_cell_serial[cell_for_i_point[i]].GetValue(element_for_i_point[i], rule_for_i_point[i]);
    }

    current_idx = 0;
    for (int i = 0; i < num_cells; i++) {
        for (int j = 0; j < num_mem_nodes_cell[i]; j++) {
            v_prev_vec[current_idx] = x_cell_serial[i][membrane_list_cell[i][j]];
            ui_prev_vec_cell[current_idx] = x_cell_serial[i][membrane_list_cell[i][j]];
            current_idx += 1;
        }
    }
    
    for (int i = 0; i < num_cells; i++) {
        v_in_cell[i] = 0.0;
        for (int j = 0; j < num_mem_nodes_cell[i]; j++) {
            v_in_cell[i][membrane_list_cell[i][j]] =  v_prev_vec[membrane_idx_for_cell[i][j]];
        }
    }

    v_in_i = 0.0;
    for (int i = 0; i < fespace_i->GetNDofs(); i++)
    {
        v_in_i[i] = v_in_cell[cell_for_i_point[i]].GetValue(element_for_i_point[i], rule_for_i_point[i]);
    }
    
        
    for (int i = 0; i < num_mem_nodes_e; i++)
    {
        v_prev_vec_e[i] = v_in_i.GetValue(membrane_elements_e[i], membrane_rule_e[i]);
        ui_prev_vec_e[i] = x_i.GetValue(membrane_elements_e[i], membrane_rule_e[i]);
    }
    ue_prev_vec_e = 0.0;
    ue_prev_vec_cell = 0.0;
    ue_new_vec_cell = 0.0;
    
    for (int i = 0; i < num_mem_nodes_total; i++)
    {
        states[i][V_idx] = v_prev_vec[i];
    }
    
    

    // Prepare for saving solution to ParaView
    std::string save_string_p = "ParaView/";
    save_string_p += simulation_name;
    
    ParaViewDataCollection pd_i("I", mesh_i);
    ParaViewDataCollection pd_m("M", mesh_i);
    ParaViewDataCollection pd_e("E", pmesh_e);
    if (save_paraview) {
        if (myid == 0) {
            pd_i.SetPrefixPath(save_string_p);
            pd_i.SetLevelsOfDetail(order);
            pd_i.SetCycle(0);
            pd_i.SetDataFormat(VTKFormat::BINARY);
            pd_i.SetHighOrderOutput(true);
            pd_i.SetTime(0.0); // set the time
            pd_i.RegisterField("$u_i$ (mV)",&x_i);
            pd_i.Save();
            
            pd_m.SetPrefixPath(save_string_p);
            pd_m.SetLevelsOfDetail(order);
            pd_m.SetCycle(0);
            pd_m.SetDataFormat(VTKFormat::BINARY);
            pd_m.SetHighOrderOutput(true);
            pd_m.SetTime(0.0); // set the time
            pd_m.RegisterField("$v$ (mV)",&v_in_i);
            pd_m.Save();
        }
    
        pd_e.SetPrefixPath(save_string_p);
        pd_e.SetLevelsOfDetail(order);
        pd_e.SetCycle(0);
        pd_e.SetDataFormat(VTKFormat::BINARY);
        pd_e.SetHighOrderOutput(true);
        pd_e.SetTime(0.0); // set the time
        pd_e.RegisterField("$u_e$ (mV)",&x_e);
        pd_e.Save();
    }

    
    
    
    // Prepare for saving V
    std::string save_v_string = "V/";
    save_v_string += simulation_name;
    save_v_string += ".txt";
    
    // Prepare for saving U
    std::string save_u_string = "U/";
    save_u_string += simulation_name;
    save_u_string += ".txt";
    
    // Save v
    ofstream wf_v;
    if (myid == 0) {
        wf_v.open(save_v_string);
        write_mean_values(wf_v, v_prev_vec, num_cells, num_mem_nodes_cell, NU_nodes, membrane_list_cell, fespace_cell_serial, mem_attr);
    }
    
    ofstream wf_u;
    if (myid == 0) {
        wf_u.open(save_u_string);
        for (int i=0; i < 9; i++) {
            wf_u << 0.0 << endl;
        }
    }
    
    
   // Set up the bilinear form for the extracellular system
   OperatorPtr A_e;
   Vector B_e, X_e;
   ParBilinearForm *a_e = new ParBilinearForm(fespace_e);
   a_e->AddDomainIntegrator(new DiffusionIntegrator(sigma_e_coeff));
   a_e->Assemble();
   b_e->Assemble();
   a_e->FormLinearSystem(ess_tdof_list_e, x_e, *b_e, A_e, X_e, B_e);

    

   // Set up the bilinear form for the intracellular systems
   OperatorPtr A_i[num_cells];
   Vector B_i[num_cells];
   Vector X_i[num_cells];
   ParBilinearForm *a_cell[num_cells];
    
   for (int i = 0; i < num_cells; i++)
   {
       a_cell[i] = new ParBilinearForm(fespace_cell[i]);
       a_cell[i]->AddDomainIntegrator(new DiffusionIntegrator(sigma_i_coeff));
       a_cell[i]->AddBoundaryIntegrator(new BoundaryMassIntegrator(Cmdt_coeff), mem_attr);
       a_cell[i]->Assemble();
       b_cell_func_serial[i] = 0.0;
       ParGridFunction b_cell_func(pmesh_cell[i], &b_cell_func_serial[i], cell_partitioning[i]);
       GridFunctionCoefficient bdr_coeff_mem(&b_cell_func);
       b_cell[i] = new LinearForm(fespace_cell[i]);
       b_cell[i]->AddBoundaryIntegrator(new BoundaryLFIntegrator(bdr_coeff_mem), mem_attr);
       b_cell[i]->Assemble();
       a_cell[i]->FormLinearSystem(ess_tdof_list_cell[i], x_cell[i], *b_cell[i], A_i[i], X_i[i], B_i[i]);
    }
    
    
    // Set up linear solver for the extracellular system
    HypreBoomerAMG M_e;
    M_e.SetAggressiveCoarsening(0);
    M_e.SetStrengthThresh(0.5);
    M_e.SetRelaxType(18);
    M_e.SetPrintLevel(0);
    CGSolver extracellularsolver(MPI_COMM_WORLD);
    extracellularsolver.SetPrintLevel(0);
    extracellularsolver.SetMaxIter(50000);
    extracellularsolver.SetRelTol(sqrt(1e-12));
    extracellularsolver.SetAbsTol(0);
    extracellularsolver.SetPreconditioner(M_e);
    extracellularsolver.SetOperator(*A_e);
    HypreBoomerAMG* M_i[num_cells];
    CGSolver* cellsolver[num_cells];
    for (int i = 0; i < num_cells; i++) {
        M_i[i] = new HypreBoomerAMG();
        M_i[i]->SetAggressiveCoarsening(0);
        M_i[i]->SetStrengthThresh(0.5);
        M_i[i]->SetRelaxType(18);
        M_i[i]->SetPrintLevel(0);
        cellsolver[i] = new CGSolver(MPI_COMM_WORLD);
        cellsolver[i]->SetPrintLevel(0);
        cellsolver[i]->SetMaxIter(50000);
        cellsolver[i]->SetRelTol(sqrt(1e-12));
        cellsolver[i]->SetAbsTol(0);
        cellsolver[i]->SetPreconditioner(*M_i[i]);
        cellsolver[i]->SetOperator(*A_i[i]);
    }

    // Timer
    auto start = chrono::steady_clock::now();
    auto timer_start = chrono::steady_clock::now();
    auto timer_end = chrono::steady_clock::now();
    
    
    double t = 0.0;
    int save_counter = 0;
    int save_inner_counter = 0;
    int print_inner_counter = 0;
    double remaining_time;

    cout << "Started simulation..." << endl;
    
    // Perform simulation
    for (int n = 0; n < Nt; n++)
    {
        
        // STEP 1: MEMBRANE MODEL ODEs
        for (int k = 0; k < K; k++)
        {
            #pragma omp parallel for
            for (int q = 0; q < num_mem_nodes_total; q++)
            {
                forward_rush_larsen(states[q], t, dt_ode, parameters_nodes[q]);
            }
            t += dt_ode;
        }
        
        // Update v_prev_vec from ODE solution
        for (int q = 0; q < num_mem_nodes_total; q++)
        {
            v_prev_vec[q] = states[q][V_idx];
        }
        
        
        // Define rhs for extracellular system
        for (int i = 0; i < num_cells; i++) {
            v_in_cell[i] = 0.0;
            for (int j = 0; j < num_mem_nodes_cell[i]; j++) {
                v_in_cell[i][membrane_list_cell[i][j]] =  v_prev_vec[membrane_idx_for_cell[i][j]];
            }
        }
        v_in_i = 0.0;
        for (int i = 0; i < fespace_i->GetNDofs(); i++)
        {
            v_in_i[i] = v_in_cell[cell_for_i_point[i]].GetValue(element_for_i_point[i], rule_for_i_point[i]);
        }
        
        for (int i = 0; i < num_mem_nodes_e; i++)
        {
             v_prev_vec_e[i] = v_in_i.GetValue(membrane_elements_e[i], membrane_rule_e[i]);
         }
        
    
        for (int k = 0; k < N_it; k++) {
            
            // STEP 2: INTRACELLULAR SYSTEMS
            for (int m = 0; m < M_it; m++)
            {
                // Update rhs of intracellular system
                for (int i = 0; i < num_cells; i++)
                {
                    b_cell[i] = new LinearForm(fespace_cell[i]);
                    
                    b_cell_func_serial[i] = 0.0;
                    for (int j = 0; j < num_mem_nodes_cell[i]; j++) {
                        b_cell_func_serial[i][membrane_list_cell[i][j]] = (Cm/dt)*(v_prev_vec[membrane_idx_for_cell[i][j]] + ue_prev_vec_cell[membrane_idx_for_cell[i][j]]);
                    }

                    b_cell_func[i] = new ParGridFunction(pmesh_cell[i], &b_cell_func_serial[i], cell_partitioning[i]);

                    GridFunctionCoefficient bdr_coeff_mem(b_cell_func[i]);
                
                    b_cell[i]->AddBoundaryIntegrator(new BoundaryLFIntegrator(bdr_coeff_mem), mem_attr);
                    b_cell[i]->Assemble();
                }

                for (int i = 0; i < num_cells; i++)
                {
                    // Form linear system
                    a_cell[i]->FormLinearSystem(ess_tdof_list_cell[i], x_cell[i], *b_cell[i], A_i[i], X_i[i], B_i[i]);
                
                    // Solve the linear system A X = B for the intracellular part.
                    cellsolver[i]->Mult(B_i[i], X_i[i]);

                    // Recover the solution as a finite element grid function.
                    a_cell[i]->RecoverFEMSolution(X_i[i], *b_cell[i], x_cell[i]);

                    // Update intracellular solution vectors
                    x_cell[i].GetTrueDofs(x_cell_vec[i]);
                    
                    delete b_cell[i];
                    delete b_cell_func[i];
                }
            }

            
            // Collect intracellular solution
            for (int i = 0; i < fespace_i->GetNDofs(); i++)
            {
                x_i[i] = x_cell[cell_for_i_point[i]].GetValue(element_for_i_point[i], rule_for_i_point[i]);
            }
                
            // STEP 3: EXTRACELLULAR SYSTEM
            for (int i = 0; i < num_mem_nodes_e; i++)
            {
                ui_prev_vec_e[i] = x_i.GetValue(membrane_elements_e[i], membrane_rule_e[i]);
            }
                
            b_e_mem_vec = ui_prev_vec_e;
            b_e_mem_vec -= v_prev_vec_e;
            b_e_mem_vec -= ue_prev_vec_e;
            b_e_mem_vec *= (Cm/dt);
            b_e_func.SetSubVector(membrane_list_e, b_e_mem_vec);
            b_e->Assemble();
                
            a_e->FormLinearSystem(ess_tdof_list_e, x_e, *b_e, A_e, X_e, B_e);
            extracellularsolver.Mult(B_e, X_e);
                
            // Recover the solution as a finite element grid function.
            a_e->RecoverFEMSolution(X_e, *b_e, x_e);
                
                
            // STEP 4: UPDATE SOLUTIONS
            current_idx = 0;
            for (int i = 0; i < num_cells; i++) {
                for (int j = 0; j < num_mem_nodes_cell[i]; j++) {
                    ue_prev_vec_cell[current_idx] = x_e.GetValue(membrane_elements_cell[i][j], membrane_rule_cell[i][j]);
                    current_idx += 1;
                }
            }
            x_e.GetSubVector(membrane_list_e, ue_prev_vec_e);
        }
        
        current_idx = 0;
        for (int i = 0; i < num_cells; i++) {
            for (int j = 0; j < num_mem_nodes_cell[i]; j++) {
                ue_new_vec_cell[current_idx] = x_e.GetValue(membrane_elements_cell[i][j], membrane_rule_cell[i][j]);
                ui_prev_vec_cell[current_idx] = x_cell[i][membrane_list_cell[i][j]];
                current_idx += 1;
            }
        }
        v_prev_vec.SetVector(ui_prev_vec_cell, 0);
        v_prev_vec -= ue_new_vec_cell;
        
        // Save V
        if (myid == 0) {
            write_mean_values(wf_v, v_prev_vec, num_cells, num_mem_nodes_cell, NU_nodes, membrane_list_cell, fespace_cell_serial, mem_attr);
        }
        
        
        
        // Save U
        if (myid == 0) {
            write_mean_values(wf_u, ue_new_vec_cell, num_cells, num_mem_nodes_cell, NU_nodes, membrane_list_cell, fespace_cell_serial, mem_attr);
        }
        
        
        // Save data in the ParaView format
        save_inner_counter += 1;
        if (save_inner_counter == save_limit)
        {
            save_inner_counter = 0;
            save_counter +=1;
            if (save_paraview) {
                if (myid == 0) {
                    pd_i.SetCycle(save_counter);
                    pd_i.SetTime((n+1)*dt);
                    pd_i.Save();
                    pd_m.SetCycle(save_counter);
                    pd_m.SetTime((n+1)*dt);
                    pd_m.Save();
                }
                pd_e.SetCycle(save_counter);
                pd_e.SetTime((n+1)*dt);
                pd_e.Save();
            }
        }
        
        if (myid == 0) {
            // Print remaining simulation time
            print_inner_counter += 1;
            if (print_inner_counter == print_limit)
            {
                print_inner_counter = 0;
                timer_end = chrono::steady_clock::now();
                remaining_time = (chrono::duration_cast<chrono::milliseconds>(timer_end-timer_start).count()/1000.0)*(Nt-n)/print_limit;
                if (remaining_time > 86400*2) {
                    cout << remaining_time/86400 << " days remaining" << endl;
                } else if (remaining_time > 3600) {
                    cout << remaining_time/3600 << " hours remaining" << endl;
                } else if (remaining_time > 60) {
                    cout << remaining_time/60 << " min remaining" << endl;
                } else {
                    cout << remaining_time << " sec remaining" << endl;
                }
                
                timer_start = chrono::steady_clock::now();
            }
        }
        

        // Update the v solution in the ODE model
        for (int q = 0; q < num_mem_nodes_total; q++)
        {
            states[q][V_idx] = v_prev_vec[q];
        }

    }

    auto end = chrono::steady_clock::now();
    
    if (myid == 0) {
        cout << "Elapsed time : " << chrono::duration_cast<chrono::milliseconds>(end-start).count()/1000.0 << " sec" << endl;
    }
    
    if (myid == 0) {
        wf_v.close();
        wf_u.close();
    }
    
    if (save_paraview) {
        if (myid == 0) {
            pd_i.SetCycle(save_counter);
            pd_i.SetTime(Tstop);
            pd_i.Save();
            pd_m.SetCycle(save_counter);
            pd_m.SetTime(Tstop);
            pd_m.Save();
        }
        pd_e.SetCycle(save_counter);
        pd_e.SetTime(Tstop);
        pd_e.Save();
        
    }
    

    
   // Free the used memory.
   delete a_e;
   delete b_e;
   delete fespace_i;
   delete fespace_e;
   delete fec_i;
   delete fec_e;
   delete mesh_i;
   delete mesh_e;
   delete[] states;


    for (int i = 0; i < num_cells; i++)
    {
        delete mesh_cell[i];
        delete a_cell[i];
        delete fespace_cell[i];
        if (order > 0) { delete fec_cell[i];}
    }

   return 0;
}





int load_number(std::string load_string)
{
    // Set up file
    load_string += ".txt";
    string line;
    ifstream in_file(load_string);
    getline(in_file, line);
    in_file.close();
    return std::stoi(line);
}

void load_number_cells(std::string load_string, int *some_number, int num_cells)
{
    // Set up file
    load_string += ".txt";
    string line;
    ifstream in_file(load_string);
    for (int i = 0; i < num_cells; i++) {
        getline(in_file, line);
        some_number[i] = std::stoi(line);
    }
    in_file.close();
}

void load_array(std::string load_string, Array<int> &some_list, int list_length)
{
    // Set up file
    load_string += ".txt";
    string line;
    ifstream in_file(load_string);
    some_list.SetSize(list_length);
    for (int i = 0; i < list_length; i++) {
        getline(in_file, line);
        some_list[i] = std::stoi(line);
    }
    in_file.close();
}

void load_array_cells(std::string load_string, Array<int> *some_list, int *list_length, int num_cells)
{
    // Set up file
    load_string += ".txt";
    string line;
    ifstream in_file(load_string);
    for (int i = 0; i < num_cells; i++) {
        some_list[i].SetSize(list_length[i]);
        for (int j = 0; j < list_length[i]; j++) {
            getline(in_file, line);
            some_list[i][j] = std::stoi(line);
        }
    }
    in_file.close();
}

void load_integration_rule(std::string load_string, IntegrationRule &some_list, int list_length)
{
    // Set up file
    std::string load_string_x = load_string;
    std::string load_string_y = load_string;
    std::string load_string_z = load_string;
    load_string_x += "_x.txt";
    load_string_y += "_y.txt";
    load_string_z += "_z.txt";

    string line;
    ifstream in_file_x(load_string_x);
    ifstream in_file_y(load_string_y);
    ifstream in_file_z(load_string_z);
    some_list.SetSize(list_length);
    for (int i = 0; i < list_length; i++) {
        //cout << i << endl;
        getline(in_file_x, line);
        some_list.IntPoint(i).x = std::stod(line);
        getline(in_file_y, line);
        some_list.IntPoint(i).y = std::stod(line);
        getline(in_file_z, line);
        some_list.IntPoint(i).z = std::stod(line);
        
    }
    in_file_x.close();
    in_file_y.close();
    in_file_z.close();
}

void load_integration_rule_cells(std::string load_string, IntegrationRule *some_list, int *list_length, int num_cells)
{
    // Set up file
    std::string load_string_x = load_string;
    std::string load_string_y = load_string;
    std::string load_string_z = load_string;
    load_string_x += "_x.txt";
    load_string_y += "_y.txt";
    load_string_z += "_z.txt";

    string line;
    ifstream in_file_x(load_string_x);
    ifstream in_file_y(load_string_y);
    ifstream in_file_z(load_string_z);
    for (int i = 0; i < num_cells; i++) {
        some_list[i].SetSize(list_length[i]);
        for (int j = 0; j < list_length[i]; j++) {
            getline(in_file_x, line);
            some_list[i].IntPoint(j).x = std::stod(line);
            getline(in_file_y, line);
            some_list[i].IntPoint(j).y = std::stod(line);
            getline(in_file_z, line);
            some_list[i].IntPoint(j).z = std::stod(line);
        }
    }
    in_file_x.close();
    in_file_y.close();
    in_file_z.close();
}



void write_mean_values(std::ofstream &wf, Vector sol_vec, int num_cells, int* num_mem_nodes_cell, Array<int> NU_nodes, Array<int> *membrane_list_cell, FiniteElementSpace **fespace_cell_serial, Array<int> mem_attr) {

    
    // Compare areas of the two cells
    GridFunction x_cell_soma[num_cells];
    GridFunction x_cell_d_small[num_cells];
    GridFunction x_cell_d_medium[num_cells];
    GridFunction x_cell_d_large[num_cells];
    GridFunction x_cell_AIS[num_cells];
    GridFunction x_cell_pAIS[num_cells];
    GridFunction x_cell_RN[num_cells];
    GridFunction x_cell_collateral[num_cells];
    GridFunction x_cell_full[num_cells];
    GridFunction f_cell_soma[num_cells];
    GridFunction f_cell_d_small[num_cells];
    GridFunction f_cell_d_medium[num_cells];
    GridFunction f_cell_d_large[num_cells];
    GridFunction f_cell_AIS[num_cells];
    GridFunction f_cell_pAIS[num_cells];
    GridFunction f_cell_RN[num_cells];
    GridFunction f_cell_collateral[num_cells];
    GridFunction f_cell_full[num_cells];
    
    ConstantCoefficient one_coeff(1.0);
    
    for (int i=0; i < num_cells; i++) {
        x_cell_soma[i].SetSpace(fespace_cell_serial[i]);
        x_cell_d_small[i].SetSpace(fespace_cell_serial[i]);
        x_cell_d_medium[i].SetSpace(fespace_cell_serial[i]);
        x_cell_d_large[i].SetSpace(fespace_cell_serial[i]);
        x_cell_AIS[i].SetSpace(fespace_cell_serial[i]);
        x_cell_pAIS[i].SetSpace(fespace_cell_serial[i]);
        x_cell_RN[i].SetSpace(fespace_cell_serial[i]);
        x_cell_collateral[i].SetSpace(fespace_cell_serial[i]);
        x_cell_full[i].SetSpace(fespace_cell_serial[i]);
        f_cell_soma[i].SetSpace(fespace_cell_serial[i]);
        f_cell_d_small[i].SetSpace(fespace_cell_serial[i]);
        f_cell_d_medium[i].SetSpace(fespace_cell_serial[i]);
        f_cell_d_large[i].SetSpace(fespace_cell_serial[i]);
        f_cell_AIS[i].SetSpace(fespace_cell_serial[i]);
        f_cell_pAIS[i].SetSpace(fespace_cell_serial[i]);
        f_cell_RN[i].SetSpace(fespace_cell_serial[i]);
        f_cell_collateral[i].SetSpace(fespace_cell_serial[i]);
        f_cell_full[i].SetSpace(fespace_cell_serial[i]);
        x_cell_full[i] = 1.0;
        x_cell_soma[i] = 0.0;
        x_cell_d_small[i] = 0.0;
        x_cell_d_medium[i] = 0.0;
        x_cell_d_large[i] = 0.0;
        x_cell_AIS[i] = 0.0;
        x_cell_pAIS[i] = 0.0;
        x_cell_RN[i] = 0.0;
        x_cell_collateral[i] = 0.0;
        f_cell_soma[i] = 0.0;
        f_cell_d_small[i] = 0.0;
        f_cell_d_medium[i] = 0.0;
        f_cell_d_large[i] = 0.0;
        f_cell_AIS[i] = 0.0;
        f_cell_pAIS[i] = 0.0;
        f_cell_RN[i] = 0.0;
        f_cell_collateral[i] = 0.0;
        f_cell_full[i] = 0.0;
    }
    
    
    double A_full[num_cells];
    double A_soma[num_cells];
    double A_d_small[num_cells];
    double A_d_medium[num_cells];
    double A_d_large[num_cells];
    double A_AIS[num_cells];
    double A_pAIS[num_cells];
    double A_RN[num_cells];
    double A_collateral[num_cells];
    
    double F_full[num_cells];
    double F_soma[num_cells];
    double F_d_small[num_cells];
    double F_d_medium[num_cells];
    double F_d_large[num_cells];
    double F_AIS[num_cells];
    double F_pAIS[num_cells];
    double F_RN[num_cells];
    double F_collateral[num_cells];
    
    
    int current_idx = 0;
    for (int i = 0; i < num_cells; i++) {
        BilinearForm *fg_m = new BilinearForm(fespace_cell_serial[i]);
        fg_m->AddBoundaryIntegrator(new BoundaryMassIntegrator(one_coeff), mem_attr);
        fg_m->Assemble();
      
        for (int j=0; j < num_mem_nodes_cell[i]; j++) {
            f_cell_full[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            if (NU_nodes[current_idx] == 4) {
                x_cell_d_large[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_d_large[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            } else if (NU_nodes[current_idx] == 5) {
                x_cell_d_medium[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_d_medium[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            } else if (NU_nodes[current_idx] == 6) {
                x_cell_d_small[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_d_small[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            } else if (NU_nodes[current_idx] == 7) {
                x_cell_soma[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_soma[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            } else if (NU_nodes[current_idx] == 8) {
                x_cell_AIS[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_AIS[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            } else if (NU_nodes[current_idx] == 9) {
                x_cell_pAIS[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_pAIS[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            } else if (NU_nodes[current_idx] == 10) {
                x_cell_RN[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_RN[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            } else if (NU_nodes[current_idx] == 11) {
                x_cell_collateral[i][membrane_list_cell[i][j]] = 1.0;
                f_cell_collateral[i][membrane_list_cell[i][j]] = sol_vec[current_idx];
            }
            current_idx++;
        }
        
        A_full[i] = (*fg_m).InnerProduct(x_cell_full[i], x_cell_full[i]);
        A_d_large[i] = (*fg_m).InnerProduct(x_cell_d_large[i], x_cell_full[i]);
        A_d_medium[i] = (*fg_m).InnerProduct(x_cell_d_medium[i], x_cell_full[i]);
        A_d_small[i] = (*fg_m).InnerProduct(x_cell_d_small[i], x_cell_full[i]);
        A_soma[i] = (*fg_m).InnerProduct(x_cell_soma[i], x_cell_full[i]);
        A_AIS[i] = (*fg_m).InnerProduct(x_cell_AIS[i], x_cell_full[i]);
        A_pAIS[i] = (*fg_m).InnerProduct(x_cell_pAIS[i], x_cell_full[i]);
        A_RN[i] = (*fg_m).InnerProduct(x_cell_RN[i], x_cell_full[i]);
        A_collateral[i] = (*fg_m).InnerProduct(x_cell_collateral[i], x_cell_full[i]);
        
        F_d_large[i] = (*fg_m).InnerProduct(f_cell_d_large[i], x_cell_full[i]);
        F_d_medium[i] = (*fg_m).InnerProduct(f_cell_d_medium[i], x_cell_full[i]);
        F_d_small[i] = (*fg_m).InnerProduct(f_cell_d_small[i], x_cell_full[i]);
        F_soma[i] = (*fg_m).InnerProduct(f_cell_soma[i], x_cell_full[i]);
        F_AIS[i] = (*fg_m).InnerProduct(f_cell_AIS[i], x_cell_full[i]);
        F_pAIS[i] = (*fg_m).InnerProduct(f_cell_pAIS[i], x_cell_full[i]);
        F_RN[i] = (*fg_m).InnerProduct(f_cell_RN[i], x_cell_full[i]);
        F_collateral[i] = (*fg_m).InnerProduct(f_cell_collateral[i], x_cell_full[i]);
        F_full[i] = (*fg_m).InnerProduct(f_cell_full[i], x_cell_full[i]);
        
        
        wf << F_full[i]/A_full[i] << endl;
        wf << F_AIS[i]/A_AIS[i] << endl;
        wf << F_soma[i]/A_soma[i] << endl;
        wf << F_d_large[i]/A_d_large[i] << endl;
        wf << F_d_medium[i]/A_d_medium[i] << endl;
        wf << F_d_small[i]/A_d_small[i] << endl;
        wf << F_pAIS[i]/A_pAIS[i] << endl;
        wf << F_RN[i]/A_RN[i] << endl;
        wf << F_collateral[i]/A_collateral[i] << endl;
        
        delete fg_m;
        
    }
    
    
}






    
