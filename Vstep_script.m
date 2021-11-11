initialise_df

%%
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

t0 = 1e-6;
tmax = 1e2;
%%
sol_stepV_int0 = stepV(soleq.ion, 0.02, 0, t0, tmax, 400);
sol_stepV_int0p001 = stepV(soleq.ion, 0.02, 0.001, t0, tmax, 400);
sol_stepV_int0p01 = stepV(soleq.ion, 0.02, 0.01, t0, tmax, 400);
sol_stepV_int0p1 = stepV(soleq.ion, 0.02, 0.1, t0, tmax, 400);
sol_stepV_int1 = stepV(soleq.ion, 0.02, 1, t0, tmax, 400);

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