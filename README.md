# IFR Successive Linearization
This repository contains work of my research assistant job at the IFR of the University of Stuttgart. The research is about successive linearization and the analysis resp. utilization
of the linearization error dynamics.

Contact: eliasniepoetter@gmail.com

## Overview
- analysis_operating_point_vdp
    - operating_point_classic_vdp.m: perform standard operating point/setpoint computation
    - operating_point_modified_vdp.m: perform operating point/setpoint computation for a modified Van der Pol
- Analysis Operating Point VdP
    - error_dynamics_linearization_point.m: mainly visualization of the error dynamics
    - linearization_error_optimimzation.m: setpoint optimization with equilibrium constraint
- experiment_van_der_pol
    - the main.m performs a full simulation with successive linearization control scheme
    - the user can select different setpoint generation methods
