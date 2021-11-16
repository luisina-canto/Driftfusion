initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%% Set up
t0 = 1e-8;                              % Initial time
tmax = 1e2;                             % Relaxation time of the current transient
tpoints = 400;                          % Time points
DeltaV = 0.2;                           % Change in voltage
intsarr = 0;     % intensity array
%% Does voltage jump
for i = 1:length(intsarr)
    sol_stepV(i) = stepV(soleq.ion, DeltaV, intsarr(i), t0, tmax, tpoints);
end
%% Analysis
for i = 1:length(intsarr)
    J = dfana.calcJ(sol_stepV(i));
    Jt = J.tot(:, 0);
    
    %VstepIsana(sol_stepV(i), 1);
    figure(51)
    hold on
    figure(52)
    hold on
    figure(53)
    hold on
end
%% Tidy up plots
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