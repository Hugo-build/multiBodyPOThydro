The repository is to set up multiple individual bodies from an entire floater of a 15MW floating wind turbine.
The .GDF file and .POT file within the "wamitInp" folder is to be used for running "Wwmit" for obtaining
the potential flow theory based radiation and diffraction library. The .5p and .pnl file from "Wamit" is mainly used.

The detailed methodology is found at 
>Ma, Y., 2022. Novel modeling and fatigue analysis for early-phase design of a 15-mw fowt (Master's thesis, NTNU). URL:https://ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/3016401


main_voturnUS.m --> is the main file for the section division and calculation

result_multi_tot.m --> is the file for comparing the different multi-body models with the rigid body model