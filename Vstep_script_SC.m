initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%% Input parameters for the current transient
t0 = 1e-10;              % Initial time after zero. The voltage steps at t = 0 and reaches V0 + DeltaV at t = t0
tmax = 10;              % End time to capture to
tpoints = 1000;         % Number of time points for the voltage step (log spacing)
DeltaV_SC = 10e-3;      % Voltage step
intsarr = [0, 1];       % Array of light intensities at which the solution will be calculated

%%
for i = 1:length(intsarr) 
    sol_SC(i) = lightonRs(soleq.ion, intsarr(i), -1, 1, 1e6, 10);       % These two lines are necessary to get close to zero current condition
    sol_SC(i) = RsToClosedCircuit(sol_Rs1e6(i));
    sol_stepV_SC(i) = stepV(sol_SC(i), DeltaV_SC, t0, tmax, tpoints);
end

%% Analysis SC
for i = 1:length(intsarr)
    VstepISana(sol_SC(i), sol_stepV_SC(i), 1);
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

%% Legend
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