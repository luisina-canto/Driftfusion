initialise_df

%%
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%%
sol_stepV_int0 = stepV(soleq.ion, 0.02, 0, 1e-3, 1e2, 400);
sol_stepV_int0p001 = stepV(soleq.ion, 0.02, 0.001, 1e-3, 1e2, 400);
sol_stepV_int0p01 = stepV(soleq.ion, 0.02, 0.01, 1e-3, 1e2, 400);
sol_stepV_int0p1 = stepV(soleq.ion, 0.02, 0.1, 1e-3, 1e2, 400);
sol_stepV_int1 = stepV(soleq.ion, 0.02, 1, 1e-3, 1e2, 400);

%%
VstepIsana(sol_stepV_int0, 1);
figure(542)
hold on
VstepIsana(sol_stepV_int0p001, 1);
VstepIsana(sol_stepV_int0p01, 1);
VstepIsana(sol_stepV_int0p1, 1);
VstepIsana(sol_stepV_int1, 1);
figure(542)
hold off