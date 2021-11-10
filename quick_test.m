input_csv = 'Input_files/spiro_mapi_tio2.csv';
par = pc(input_csv);
soleq = equilibrate(par);
sol_1sun = changeLight(soleq.ion, 1, 100, 1);
[sol_1sun_Rs1e6, Voc] = findVocDirect(sol_1sun, 1, true, 20);

%% Switch off series resistance and apply Vapp = Voc
% Rather than sweeping forward from equilibrium here I have used the
% solution from findVocDirect, switched off the series resistance and
% applied Vapp = Voc. Since this is something we may wish to do often I
% have 
sol_closed = RsToClosedCircuit(sol_1sun_Rs1e6);

%% Plot
dfplot.ELx(sol_closed)