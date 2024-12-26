# teclistamab
Code for "Using Virtual Patients to Evaluate Dosing Strategies for Bispecifics with a Bell-Shaped Efficacy Curve"
1. pop_PK_95percConfidence.m
   - Generates green curves in Figure 3 when run at dose of 0.72
   - Generates red curves in Figure 3 when run at dose of 1.5
   - Generates Figure 7A and 7B when run at dose of 0.72
   - Generates Figures 7C and 7D when run at dose of 1.5
   - Note that the code requires VPs_dose0pt72.mat to run at dose of 0.72, and VPs_dose1pt5.mat to run at dose of 1.5
   - However, both of these files were too large to store on GitHub. They are available for download on Google Drive:
     Dose of 0.72: https://drive.google.com/file/d/17-1_qh7lbqopBSgEuBYz-CVsS7o4UyXp/view?usp=drive_link 
     Dose of 1.5:  https://drive.google.com/file/d/1x4mp-Aprx_TO9jYKuP15dS6jELiHRrvo/view?usp=drive_link 

2. trimer_mm_sensitivity_final.m
   - Generates the data in Figure 4's sensitivty analysis

3. sweep_quartiles_5_params_final.m
   - Generates the data in Figure 5-6 and 8-10
