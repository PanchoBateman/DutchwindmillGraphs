# ComputeDWMGraphsJCF
Matlab function which computes the Jordan Canonical Form and the transformation matrices for any Directed Dutch windmill graph of eta cycles of kappa nodes plus a common vertex.
Syntax : [A,L,J,P,eta, kappa,nb_jordan,sizejordan,coljordan]=ComputeDWMGraphsJCF
Type help ComputeDWMGraphsJCF to obtain the description of the variables

# SecondOrderConsensusStateMatrixJCF.m
Matlab file (to run) which computes and displays the Jordan Canonical form of the second-order consensus state space matrix 
for a precomputed Laplacian matrix. If desired, this program also computes the Laplacian matrix of a Directed Dutchwindmill graph
and the associated second-order consensus state matrix. 

# Dutchwindmillxy.mat
Data files containing the adjacency, Laplacian, transformation matrix, JCF, etc. for a Directed Dutch Windmill graph with x cycles and y vertices + 1 common vertex for each cycle.
This matfiles have been obtained with the ComputeDWMGraphsJCF.m file and can be used with the SecondOrderConsensusStateMatrixJCF program.

# Sim_Second_Order_Consensus.slx 
Simulink model which provides the time responses for a second-order consensus. It requires the matrices obtained with the SecondOrderConsensusStateMatrixJCF.m program.

  

