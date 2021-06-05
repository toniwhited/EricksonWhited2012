The main executable file is CrSc_mod.m.      You can change this to change the seeds. Search for the line "seed = 1;"
The main estimation file is genBetaWD_cr.m   To change between the GEARY and other starting values algorithm, go to line 640 and toggle the flag between 0 and 1.

To replicate Table 6, you will have to run CrSc_mod.m two times. For both runs, let the seed go from 1 to 2. For the first one, use the GEARY starting value. For the other use the GEARY/OLS starting values.
