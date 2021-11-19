initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2_IS.csv');
soleq = equilibrate(par);

%% Input parameters for the current transient
t0 = 1e-9;
tmax = 1e2;
tpoints = 600;
DeltaV_SC = 20e-3;
intsarr = [0, 1];

%%
for i = 1:length(intsarr)
    if intsarr(i) == 0
        sol_SC(i) = soleq.ion;
    else      
        sol_SC(i) = lightonRs(soleq.ion, intsarr(i), -100, 1, 0, 10);
    end
    sol_stepV_SC(i) = stepV(sol_SC(i), DeltaV_SC, t0, tmax, tpoints);
end

%% Analysis SC
for i = 1:length(intsarr)
    VstepISana(sol_stepV_SC(i), 1);
    figure(51)
    hold on
    figure(52)
    hold on
    figure(53)
    hold on
    figure(54)
    hold on
    figure(55)
    hold on
end

%% Legends
legendCell = cellstr(num2str(intsarr', '%-d Sun'));
figure(51)
legend(legendCell)
hold off
figure(52)
legend(legendCell)
hold off
figure(53)
legend(legendCell)
hold off
figure(54)
legend(legendCell)
hold off
figure(55)
legend(legendCell)
hold off