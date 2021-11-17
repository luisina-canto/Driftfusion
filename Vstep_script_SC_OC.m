initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%% Input parameters for the current transient
t0 = 1e-12;
tmax = 1e2;
tpoints = 800;
DeltaV_SC = 0.05;
DeltaV_OC = 0.05;
intsarr = [1];
%%
for i = 1:length(intsarr)
    if intsarr(i) == 0
        sol_SC(i) = soleq.ion;
    else      
        sol_SC(i) = lightonRs(soleq.ion, intsarr(i), -100, 1, 0, 10);
    end
    sol_stepV_SC(i) = stepV(sol_SC(i), DeltaV_SC, t0, tmax, tpoints);
end

%% Obtain open circuit solutions
for i = 1:length(intsarr)
    sol_Rs1e6(i) = lightonRs(soleq.ion, intsarr(i), -100, 1, 1e6, 100);
    sol_OC(i) = RsToClosedCircuit(sol_Rs1e6(i));
    sol_stepV_OC(i) = stepV(sol_OC(i), DeltaV_OC, t0, tmax, tpoints);
end

%% Analysis SC
for i = 1:length(intsarr)
    VstepISana(sol_stepV_SC(i), 1);
    figure(51)
    hold on
    figure(52)
    hold onw
    figure(53)
    hold on
end

%% Analysis OC
for i = 1:length(intsarr)
    VstepISana(sol_stepV_OC(i), 1);
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
xlim([t0, tmax])
hold off
figure(52)
legend(legendCell)
xlim([1/tmax, 1/t0])
hold off
figure(53)
legend(legendCell)
xlim([1/tmax, 1/t0])
hold off