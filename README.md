# teclistamab
Code for "Using Virtual Patients to Evaluate Dosing Strategies for Bispecifics with a Bell-Shaped Efficacy Curve"
1. trimer_teclistimab._vs_data.m
   - Generates fits to data shown in Figure 2

2. individual_VP_trajectories_noDisplay.m
   - Generates individual trajectories in Figure 4, S4 and S5 when run at dose of 2
   - Reads in data from output_teclistimab_VCT_parallel.mat

3. trimer_mm_sensitivity_final.m
   - Generates the data in Figure 3's sensitivity analysis

4. protocol_sweep_VPs_subgroup.m
   - Generates the data in Figure 5-6 and S1-S4
   - The output of running this code is found in output_teclistimab_VCT_parallel.mat
