initialise_df

%% Load parameters
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);

%% Input parameters for the current transient
t0 = 1e-11;
tmax = 10;
tpoints = 1000;
DeltaV_OC = 20e-3;
intsarr = [0, 1];

%% Obtain open circuit solutions
for i = 1:length(intsarr)
    if intsarr(i) == 0
        sol_OC(i) = soleq.ion;
    else
        sol_Rs1e6(i) = lightonRs(soleq.ion, intsarr(i), -100, 1, 1e6, 100);
        sol_OC(i) = RsToClosedCircuit(sol_Rs1e6(i));
    end
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