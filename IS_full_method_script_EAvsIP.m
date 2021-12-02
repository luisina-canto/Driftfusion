% UNIT_TEST - Runs various simulations in order to test the functioning of Driftfusion and many other implemented functions
% For profiling the time spent running each function run 'profile on' before running the tests and
% 'profile viewer' after.
% Using Matlab's Coverage Reports, the obsolete and unused code can be easily spotted.
%
% Syntax:  runtests('unit_test')
%
% Other m-files required: pc, equilibrate, doJV, explore
% Subfunctions: none
% MAT-files required: none
%
% See also df.
%
% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% prepare solutions and indirectly test equilibrate and genIntStructs
% functions
initialise_df

input_csv = 'Input_files/spiro_mapi_tio2.csv';

% par = pc(varargin)
par = pc(input_csv);

% Other parameters to change:
% HTL doping level: par.EF0(1)      % (eV)
% For Ohmic contact- also change Phi_left fot same value
% ETL doping level: par.EF0(end)    % (eV)
% For Ohmic contact- also change Phi_right fot same value
% Active layer thickness: par.d(3)  % (cm)
% Left-hand workfunction: par.Phi_left      % (eV)
% Right-hand workfunction: par.Phi_right    % (eV)

EF0_arr = [0,1,2,3,4];
EA_arr = [-2.9,-2.95,-3,-3.05,-3.1];
IP_arr = [-5.0,-5.05,-5.1,-5.15,-5.2];

%-2.9,-2.95,-3,-3.05,-3.1
%-5.0,-5.05,-5.1,-5.15,-5.2

% Ncat_arr = [1e15, 1e16, 1e17, 1e18, 1e19];
%EF0_arr = -4.7;
Nani_arr = [1e16,1e17];

% Calculating tghe mobility which varies with the reciprocal of ion
% concentration to give a constant conductivity
%mu_carr = 1e19*1e-10./Ncat_arr;

for i = 1:length(EF0_arr)
    par_struct(i) = par;
    
    %par_struct(i).mu_c(3) = mu_carr(i);
    %par_struct(i).Ncat = [Nani_arr(i), Nani_arr(i), Nani_arr(i), Nani_arr(i), Nani_arr(i)];
    
    par_struct(i).Phi_IP(1) = IP_arr(i);
    par_struct(i).Phi_EA(1) = EA_arr(i);
    %par_struct(i).Phi_left = EF0_arr(i);
    
    % Everytime you change your parameters in a script use this function:
    par_struct(i) = refresh_device(par_struct(i));
    % soleq = equilibrate(varargin)
    soleq(i) = equilibrate(par_struct(i));
end

for i = 1:length(EF0_arr)
    %% Obtain open circuit initial condition
    % sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
    sol_Rs1e6(i) = lightonRs(soleq(i).ion, 0.5, -100, 1, 1e6, 100);
    sol_OC(i) = RsToClosedCircuit(sol_Rs1e6(i));
    'This_works'
end

%% Scripts IS, helper and analysis

startFreq = 1e7;  % Maybe change to 1e7?
endFreq = 1e-2;
Freq_points = 25;
deltaV = 2e-3;
frozen_ions = false;
demodulation = true;
do_graphics = false;
for i = 1:length(EF0_arr)
    % IS_results = IS_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics)
    sol_IS_SC(i) = IS_script(soleq(i).ion, startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics);
    
    % IS_results = IS_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics)
    sol_IS_OC(i) = IS_script(sol_OC(i), startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics);
end

%% Analysis plots
% IS_script_exporter(prefix, IS_results)
ea_ip = [EA_arr;IP_arr];
for i = 1:length(EF0_arr)
    %IS_script_exporter('unit_testing_deleteme', sol_IS_OC(i))

    % IS_script_plot_impedance(IS_results)
    IS_script_plot_impedance_2(sol_IS_OC(i),i,ea_ip(:,i),length(EF0_arr))
    
    % IS_script_plot_nyquist(IS_results)
    IS_script_plot_nyquist_2(sol_IS_OC(i),i,ea_ip(:,i),length(EF0_arr))

    % IS_script_plot_phase(IS_results)
    IS_script_plot_phase_2(sol_IS_OC(i),i,ea_ip(:,i),length(EF0_arr))
end 
%% Scripts IS non parallel and analysis

% % IS_results = IS_script_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, demodulation, do_graphics, save_solutions)
% % one input
% IS_dark1 = IS_script_nonparallel(soleq.ion, 1e9, 1e-2, 3, 2e-3, false, true, true, true, true);
% IS_dark2 = IS_script_nonparallel(soleq.ion, 1e9, 1e-2, 3, 2e-3, true, false, false, false, false);
% % many inputs
% IS_script_nonparallel(structs_sc, 1e9, 1e-2, 2, 2e-3, false, false, true, true, false);
% 
% % IS_list_plot(type, dir_file_name, varargin)
% IS_list_plot('capacitance', 'unit_testing_deleteme',...
%     IS_dark1, 'frozen ions',{'--r'},...
%     IS_dark2, 'sequential',{':k','LineWidth',3});
% IS_list_plot('impedance_re', 'unit_testing_deleteme',...
%     IS_dark1, 'frozen ions',{'--r'},...
%     IS_dark2, 'sequential',{':k','LineWidth',3});
% IS_list_plot('impedance_im', 'unit_testing_deleteme',...
%     IS_dark1, 'frozen ions',{'--r'},...
%     IS_dark2, 'sequential',{':k','LineWidth',3});
% IS_list_plot('impedance_abs', 'unit_testing_deleteme',...
%     IS_dark1, 'frozen ions',{'--r'},...
%     IS_dark2, 'sequential',{':k','LineWidth',3});
% IS_list_plot('phase', 'unit_testing_deleteme',...
%     IS_dark1, 'frozen ions',{'--r'},...
%     IS_dark2, 'sequential',{':k','LineWidth',3});
% IS_list_plot('nyquist', 'unit_testing_deleteme',...
%     IS_dark1, 'frozen ions',{'--r'},...
%     IS_dark2, 'sequential',{':k','LineWidth',3});

% %% Scripts EA, helper and analysis
% % 
% % EA_results = EA_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, do_graphics)
% % one input
% EA_dark1 = EA_script(soleq.ion, 1e9, 1e-2, 3, 1e-3, true, false);
% EA_dark2 = EA_script(soleq.ion, 1e9, 1e-2, 3, 1e-3, false, true);
% % many inputs
% EA_sc = EA_script(structs_sc, 1e9, 1e-2, 2, 1e-3, false, true);
% 
% % EA_script_exporter(prefix, EA_results)
% EA_script_exporter('unit_testing_deleteme', EA_sc);
% 
% % EA_list_plot(type, dir_file_name, varargin)
% EA_list_plot('2h', 'unit_testing_deleteme',...
%     EA_dark1, 'frozen_ions',{':r','LineWidth',3},...
%     EA_dark2, 'mobile ions',{'-k'});
% EA_list_plot('1h', 'unit_testing_deleteme',...
%     EA_dark1, 'frozen_ions',{':r','LineWidth',3},...
%     EA_dark2, 'mobile ions',{'-k'});
% EA_list_plot('phase', 'unit_testing_deleteme',...
%     EA_dark1, 'frozen_ions',{':r','LineWidth',3},...
%     EA_dark2, 'mobile ions',{'-k'});
% 
% % EA_script_plot_Efield(EA_results, savefig_dir)
% EA_script_plot_Efield(EA_sc, "unit_testing_deleteme")
% 
% % EA_script_plot_phase(EA_results, savefig_dir)
% EA_script_plot_phase(EA_sc, "unit_testing_deleteme");
% 
% %% Scripts SDP and helper
% 
% % sdpsol = SDP_script(sol_ini, tdwell_arr, Vjump, bias_source, bias_int, pulse_source, pulse_int, pulse_tmax, pulse_mobile_ions, solver_split_pulse)
% sdpsol1 = SDP_script(soleq.ion, logspace(-8,3,3), 0.6, 1, 0.1, 2, 5.12, 1e-3, true, true);
% sdpsol2 = SDP_script(soleq.ion, logspace(-8,3,3), 0.6, 1, 0.1, 2, 5.12, 1e-3, false, false);
% 
% % SDP_script_exporter(prefix, varargin)
% SDP_script_exporter('unit_testing_deleteme', sdpsol1, sdpsol2);
% 
% % SDP_script_plot(Jtr_time, dir_file_name, varargin)
% SDP_script_plot(1e-7, 'unit_testing_deleteme',...
%     sdpsol1, 'pulse mobile ions', {':r','LineWidth',3},...
%     sdpsol2, 'pulse frozen ions', {'-k'})
% 
% %% Protocols doIS_EA, helper and analysis
% 
% % struct_IS = doIS_EA(struct_Int, deltaV, freq, periods, tpoints_per_period, stability_timefraction, RelTol)
% is_ea_10mHz_100mV = doIS_EA(soleq.ion, 0.1, 1e-3, 20, 40, 0.5, 1e-6);
% 
% % IS_EA_struct_exporter(prefix, struct)
% IS_EA_struct_exporter('unit_testing_deleteme', is_ea_10mHz_100mV);
% 
% % coeff = EA_ana_plot(struct_IS_EA, do_graphics, local_field, demodulation, savefig_dir)
% EA_ana_plot(is_ea_10mHz_100mV, true, true, true, "unit_testing_deleteme");
% EA_ana_plot(is_ea_10mHz_100mV, false, false, false, missing);
% 
% % coeff = IS_ana_plot(s, do_graphics, demodulation)
% IS_ana_plot(is_ea_10mHz_100mV, true, true);
% IS_ana_plot(is_ea_10mHz_100mV, false, false);
% 
% J = dfana.calcJ(is_ea_10mHz_100mV).tot(:,end);
% % coeff = IS_EA_ana_demodulation(t, y, fun_type, freq)
% IS_EA_ana_demodulation(is_ea_10mHz_100mV.t', J, 'sin', is_ea_10mHz_100mV.par.V_fun_arg(3));
% 
% % coeff = IS_EA_ana_fit(t, y, fun_type, freq)
% IS_EA_ana_fit(is_ea_10mHz_100mV.t', J, 'sin', is_ea_10mHz_100mV.par.V_fun_arg(3))
% 
% % [subtracting_q_t, subtracting_q_intr_t] = IS_ana_subtracting(struct_IS)
% IS_ana_subtracting(is_ea_10mHz_100mV);
% 
% %% Protocols doLightPulse
% 
% % [sol_pulse] = doLightPulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
% doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, true);
% doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, false);
% 
% %% Protocols doSDP and Analysis anasdp
% 
% % sdpsol = doSDP(sol_ini, tdwell_arr, Vjump, pulse_int, pulse_tmax, duty, tpoints, scalefactor)
% sdpsol = doSDP(soleq.ion, [1e-3,1], 0.9, 1, 1, 5, 50, 1);
% 
% % anasdp(sdpsol, Jtr_time)
% anasdp(sdpsol, 1);
% 
% %% Protocols doSDP_alt
% 
% % sdp = doSDP_alt(sol_ini, sol_jump_is_given, tdwell_arr, Vjump, bias_source, bias_int, pulse_source, pulse_int, pulse_tmax, pulse_mobile_ions, solver_split_pulse)
% sdpsol = SDP_script(soleq.ion, logspace(-8,3,3), 0.6, 1, 0.1, 2, 5.12, 1e-3, true, true);
% 
% % anasdp(sdpsol, Jtr_time)
% anasdp(sdpsol, 1e-4);
% 
% %% Protocols doSPV and Analysis spvana
% 
% % spvsol = doSPV(sol_ini, Int, mobseti, tpoints, tmax, Rs, stabilise)
% spvsol = doSPV(soleq.ion, 1, true, 20, 1, 1e3, true);
% 
% % spvdat = spvana(spvsol)
% spvana(spvsol);

%------------- END OF CODE --------------

