par = pc('Input_files/spiro_mapi_tio2.csv');
eqm = equilibrate(par);
CV_solution = doCV(eqm.ion, 1, 0, 1.1, 0, 10e-3, 1, 211);

%% Plots
Jcomp = Current_Contributions(CV_solution);
stats = CVstats(CV_solution);