# C_code_for_numerics_in_some_quantum_phase_transitions

This repository explores some finite size scaling properties of location of turning points and peak heights of a susceptibility of entanglement entropy in transverse field Ising model and XY model(to be added) in 1D, using C. Also precision checks between computations in double and long double type are present. 



-----------------------------


` gcc -o tfim_sigma _z __ tfim_sigma _z __.c-lm`

`./tfim_sigma_z __ <h_init> <h_final> <samplesize> <N>`


-----------------------------


`gcc -o tfim_suscept_ent_ent__ tfim_suscept_ent_ent__.c -lm`

`./tfim_suscept_ent_ent__ <h_init> <h_final> <samplesize> <N>`


-----------------------------


`gcc -o tfim_suscept_ent_ent_scaling tfim_suscept_ent_ent_peak_scaling__.c -lm`

`./tfim_suscept_ent_ent_scaling`


-----------------------------


`gcc -o tfim_suscept_peak_height_scaling__ tfim_suscept_ent_ent_peak_height_scaling.c -lm`

`./tfim_suscept_peak_height_scaling__`


-----------------------------


For generating the data in `long double` computation

`gcc -o dataCompTFIM_in_LongDouble tfim_suscept_ent_ent__dataComparison.c -lm`

`./dataCompTFIM_in_LongDouble <N>`   


----------------------------



For generating the data in `double` computation

`gcc -o dataCompTFIM_in_Double tfim_suscept_ent_ent__inDouble_dataComparison.c -lm`

`./dataCompTFIM_in_Double <N>`


----------------------------
