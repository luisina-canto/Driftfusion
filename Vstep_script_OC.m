initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%% Input parameters for the current transient
t0 = 1e-9;              % Initial time after zero. The voltage steps at t = 0 and reaches V0 + DeltaV at t = t0
tmax = 10;              % End time to capture to
tpoints = 800;         % Number of time points for the voltage step (log spacing)
DeltaV_OC = 10e-3;      % Voltage step
intsarr = [0, 1e-3, 1e-2, 1e-1, 1];       % Array of light intensities at which the solution will be calculated

%% Obtain open circuit solutions
for i = 1:length(intsarr)
    sol_Rs1e6(i) = lightonRs(soleq.ion, intsarr(i), -100, 1, 1e6, 100);     % This step is required for dark conditions also to get close to a zero current condition.
    sol_OC(i) = RsToClosedCircuit(sol_Rs1e6(i));
    sol_stepV_OC(i) = stepV(sol_OC(i), DeltaV_OC, t0, tmax, tpoints);
end

%% Analysis OC
for i = 1:length(intsarr)
    VstepISana(sol_OC(i), sol_stepV_OC(i), 1);
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