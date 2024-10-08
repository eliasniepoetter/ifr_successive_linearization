# IFR Successive Linearization
This repository contains work of my research assistant job at the IFR of the University of Stuttgart. The research is 
about successive linearization and the analysis resp. utilization of the linearization error dynamics. For more details 
about the basic setup of the work in this repository, it is recommended to read the notes in the material folder.

Contact: eliasniepoetter@gmail.com

## Overview
**analysis_operating_point_vdp**:
- `operating_point_classic_vdp.m`: perform standard operating point/setpoint computation
- `operating_point_modified_vdp.m`: perform operating point/setpoint computation for a modified Van der Pol

**analysis_error_dynamics**:
- `error_dynamics_linearization_point.m`: mainly visualization of the error dynamics
- `linearization_error_optimimzation.m`: setpoint optimization with equilibrium constraint

**experiment_van_der_pol**:
- the `main.m` performs a full simulation with successive linearization control scheme
- the user can select different setpoint generation methods

## Analysis of the Operating Point of Van der Pol
These scripts demonstrate the classic design of a controller around an operating point.
The operating point is in general not an equilibrium. Therefore the equilibrium equations are
solved for a steady state input $u^{\star}$, which enables the shift of the equilibrium to (any) selected
operating point. For the "classic" Van der Pol with only one input, the operating point
is tight to the $x_1$ axis because full state tracking would require as many inputs as states.
The "modified" Van der Pol in the example presendet has an additional input for the ODE of 
$x_1$, namely $\dot{x_1} = x_2 + u_1$. This modification allows it to move the operating point
resp. the equilibrium to any point in the state space. However, one drawback is now that the equilibrium
states for $x_1$ and $x_2$ are coupled. Furthermore, if one takes the steady state input $u_{2}^{\star}$ as 
a parameter, the equilibrium undergoes a bifurcation.

## Analysis of the linearization error dynamics
The error_dyanmics_linearization_point.m script visualizes the error dynamics $e(x,u)$ for the Van der Pol.
The analytically computed results of the gradient of the linearization error are visualized as well.
On this basis, the linearization_error_optimization.m script implements an optimization routine
which searches for a setpoint with minimal linearization error (evaluated at the current state) while 
simultaneously ensuring an equilibrium at the new setpoint through a nonlinear constraint.

## Numerical Experiment: Successive Control of the Van der Pol 
The experiment performs a full simulation of the "classic" Van der Pol, controlled with a successive linearization and control
design. In this example, the operating point computation is not done fully correctly, which means that the affine dynamics
are simply dropped for simplicity. The user can select two different strategies for the setpoint generation.
- v1 and v2: no valid strategies to determine the setpoints
- v3: sampling-based approach
- v4: optimization based approach










