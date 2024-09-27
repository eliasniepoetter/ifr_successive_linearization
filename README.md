# IFR Successive Linearization
This repository contains work of my research assistant job at the IFR of the University of Stuttgart. The research is about successive linearization and the analysis resp. utilization
of the linearization error dynamics.

Contact: eliasniepoetter@gmail.com

## Overview
**analysis_operating_point_vdp**:
- operating_point_classic_vdp.m: perform standard operating point/setpoint computation
- operating_point_modified_vdp.m: perform operating point/setpoint computation for a modified Van der Pol

**analysis_error_dynamics**:
- error_dynamics_linearization_point.m: mainly visualization of the error dynamics
- linearization_error_optimimzation.m: setpoint optimization with equilibrium constraint

**experiment_van_der_pol**:
- the main.m performs a full simulation with successive linearization control scheme
- the user can select different setpoint generation methods

## Analysis of the Operating Point of Van der Pol
These scripts demonstrate the classic design of a controller around an operating point.
The operating point is in general not an equilibrium. Therefore the equilirbium equations are
solved for a steady state input $u^*$, which enables the shift of the equlibrium to (any) selected
operating point. For the "classic" Van der Pol with only one input, the operating point
is tight to the $x_1$ axis because full state tracking would require as many inputs as states.
The "modified" Van der Pol in the example presendet has an additional input for the ODE of 
$x_1$, namely $\dot{x_1} = x_2 + u_1$.
