function sol_stepV = stepV(sol_in, DeltaV, int, t0, tmax, tpoints)
% Alternative to jumptoV using Ilario's excellent 'sweepAndStill' function
% N is number of periods over which FT is calculated
if int ~= 0
    sol = lightonRs(sol_in, int, -1, 1, 0, 40);
else
    sol = sol_in;
end

par = sol.par;
par.tmesh_type = 'log10';
par.tmax = tmax;
par.t0 = t0;
par.tpoints = tpoints;
t = meshgen_t(par);

par.V_fun_type = 'sweepAndStill';
% COEFF = [A_start, A_end, sweep_duration]
V0 = getVend(sol_in);
par.V_fun_arg(1) = V0;
par.V_fun_arg(2) = V0 + DeltaV;
par.V_fun_arg(3) = t(2);            % Sweeps between 0 and first non-zero time point

disp('Performing voltage step')
sol_stepV = df(sol, par);
disp('complete')

end