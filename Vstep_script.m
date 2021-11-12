initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%% Set up 1
t0 = 1e-8;
tmax = 1e2;
tpoints = 400;
DeltaV = 0.2;
intsarr = [0, 1e-2, 1];
%%
for i = 1:length(intsarr)
    sol_stepV(i) = stepV(soleq.ion, DeltaV, intsarr(i), t0, tmax, tpoints);
end
%% Analysis
for i = 1:length(intsarr)
    VstepIsana(sol_stepV(i), 1);
    figure(51)
    hold on
    figure(52)
    hold on
    figure(53)
    hold on
end
%%
legendCell = cellstr(num2str(intsarr', '%-d Sun'));
figure(51)
legend(legendCell)
xlim([1e-2, 1e6])
hold off
figure(52)
legend(legendCell)
xlim([1e-2, 1e6])
hold off
figure(53)
legend(legendCell)
xlim([1e-2, 1e6])
hold off