function sol_stepV = stepV(sol_in, DeltaV, t0, tmax, tpoints)
% Alternative to jumptoV using Ilario's excellent 'sweepAndStill' function
% N is number of periods over which FT is calculated
% DELTAV = Change in voltage from the initial voltage in SOL_IN
% INT = Light intensity
% T0 = first time point after zero and time at which the voltage has changed
% TMAX = Maximum time
% TPOINTS = Number of time points
tic

par = sol_in.par;
par.tmesh_type = 'log10';
par.tmax = tmax;
par.t0 = t0;
par.tpoints = tpoints;
t = meshgen_t(par);
% Obtain maximum time step
tstepmax = max(diff(t));
par.MaxStepFactor = 0.1*tstepmax/par.tmax;

par.V_fun_type = 'sweepAndStill';
% COEFF = [A_start, A_end, sweep_duration]
V0 = getVend(sol_in);
par.V_fun_arg(1) = V0;
par.V_fun_arg(2) = V0 + DeltaV;
par.V_fun_arg(3) = t(2);            % Sweeps between 0 and first non-zero time point

disp('Performing voltage step')
sol_stepV = df(sol_in, par);
disp('complete')
toc
end