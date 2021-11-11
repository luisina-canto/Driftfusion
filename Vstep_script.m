initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%% Set up 1
t0 = 1e-6;
tmax = 1e2;
DeltaV = 0.05;
intsarr = [0, 1e-3, 1e-2, 1e-1, 1];
%%
for i = 1:length(intsarr)
    sol_stepV(i) = stepV(soleq.ion, DeltaV, intsarr(i), t0, tmax, 1000);
end
%% Analysis
for i = 1:length(intsarr)
    VstepIsana(sol_stepV(i), 1);
    figure(542)
    hold on
end
figure(542)
hold off