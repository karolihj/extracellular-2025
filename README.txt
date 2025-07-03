Code for running the simulations in the manuscript KH Jaeger and A Tveito, 'Sometimes extracellular recordings fail for good reasons', 2025.

The NeuronCable directory contains MATLAB code for running cable equation simulations of a simplified neuron. The run_simulation.m script can be used to run a simulation.

The SAN directory contains MATLAB code for running bidomain model simulations of a sinoatrial node tissue sample. The run_simulation.m script can be used to run a simulation.

The Heart-on-a-chip directory contains MATLAB code for running Kirchhoff's network model (KNM) simulations of a collection of hiPSC-CMs. The run_simulation.m script can be used to run a simulation.

The PancreaticIslet directory contains MATLAB code for running KNM simulations of an islet of pancreatic beta-cells. The run_simulation.m script can be used to run a simulation.

The Purkinje directory contains C++ code for running EMI model simulations of one or two cerebellar purkinje neurons. Simulations of one cell can be performed using EMI_one_cell.cpp and simulations of two cells can be performed using EMI_two_cells.cpp. This code uses the MFEM finite element library available at https://mfem.org.



