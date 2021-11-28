
%IS_SCRIPT_PLOT_IMPEDANCE - Represents Bode plots of impedance and capacitance from impedance spectroscopy
% Plot the value of absolute impedance magnitude |Z|, resistance Z'
% (real component), reactance Z'' (imaginary component), and apparent
% capacitance defined as \omega^-1 \Im{Z^-1} at various frequencies and
% various illumination or background voltage bias conditions.
% As a reference, the following currents are also used for plotting the
% mentioned quantities: ionic displacement current, recombination current,
% and total charge time derivative as obtained from
% IS_ana_subtracting
%
% Syntax:  IS_script_plot_impedance(IS_results)
%
% Inputs:
%   IS_RESULTS - a struct containing the most important results of the ISstep or the IS simulation
%
% Example:
%   IS_script_plot_impedance(IS_oc)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also IS_script, IS_script_ana_phase, IS_script_ana_nyquist, IS_list_plot.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% check which was the variable being explored
if numel(unique(sol_CV.Int)) > 1
    legend_text = IS_results.Int;
    legend_append = ' sun';
else
    legend_text = IS_results.Vdc;
    legend_append = ' Vdc';
end

% create a color array with one color more than necessary
jet_matrix = jet(length(legend_text) + 1);
% find the yellow (which in RGB code is 1,1,0) and remove it from the
% colors list
jet_yellow_logical = ismember(jet_matrix, [1, 1, 0], 'rows');
jet_no_yellow = jet_matrix(~jet_yellow_logical, :);
jet_no_yellow_flip = flipud(jet_no_yellow);
colormap cool;
Int_colors = colormap;

% round to two significant digits
legend_flip = round(legend_text, 2, 'significant');
% flip array and matrixes for starting from dark
legend_flip = string(flipud(legend_flip));
% add sun to numbers in legend
legend_flip = strcat(legend_flip, legend_append);
% replace zero in legend with dark
legend_flip(legend_flip=="0 sun") = "dark";

% preallocate figures handles
h = zeros(length(legend_text), 1);

% in case just a single frequency was simulated, min and max are
% coincident, so they are modified by a small percentage
xlim_array = [min(min(IS_results.Freq))*0.99, max(max(IS_results.Freq)*1.01)];

%% Output current and voltage for first light intensity (usually 0)
Vapp = dfana.calcVapp(sol_CV);  % Applied voltage
Jstruct = dfana.calcJ(sol_CV);
J = Jstruct.tot(:,1);                   % Total current at left-hand boundary


%% Plot energy level diagram at applied bias Vapp for first light intensity
Vplot = 0;
% Get corresponding time, TPLOT for VPLOT']
    
if Vplot >= 0
    tplot = Vplot/scan_rate;
elseif Vplot < 0
    tplot = ((2*Vmax)+abs(Vplot))/scan_rate; 
end
% PLot the energy level diagram at time TPLOT
%dfplot.ELx(sol_CV, tplot);
% PLot the carrier densities at time TPLOT
%dfplot.npx(sol_CV, tplot);


%% plot apparent capacitance vs frequency


leg = join(["Ncat =", string(change)]);
col = (floor(length(Int_colors)/6))*loop;

figure('Name', 'Energy Diagrams', 'NumberTitle', 'on');
    hold on
    for i = 1:length(legend_text)
        legend_flip;
        figure(2)
        
        h(i) = dfplot.ELx(sol_CV(i), tplot)
            'Color', Int_colors(col,:), 'MarkerEdgeColor', Int_colors(col, :),...
            'MarkerFaceColor', Int_colors(col, :), 'Marker', 's',...
  %          'MarkerSize',3, 'LineWidth', 1.3,'DisplayName',leg) ;
        legend;
        hold on
        % for plotting the negative capacitance
%         plot(IS_results.Freq(i, :), -IS_results.cap(i, :)',...
%             'Color', Int_colors(loop*25, :), 'MarkerEdgeColor', Int_colors(loop*25, :),...
%             'MarkerFaceColor', 'White', 'Marker', 's',...
%             'MarkerSize', 7, 'LineWidth', 1.3);
        if isfield(IS_results, 'cap_ion_disp')
            % capacitance due to ions
            %plot(IS_results.Freq(i, :), IS_results.cap_ion_disp(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
            %plot(IS_results.Freq(i, :), IS_results.cap_cat_disp(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1, 'Marker', '+');
            %plot(IS_results.Freq(i, :), IS_results.cap_ani_disp(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1, 'Marker', 'x');
            % capacitance calculated from the recombination amount
            %plot(IS_results.Freq(i, :), IS_results.cap_r(i, :)', 'Color', Int_colors(i, :), 'LineStyle', ':', 'LineWidth', 1.5);
        end
        if isfield(IS_results, 'cap_np_dt') % not present in ISstep simulation
            %plot(IS_results.Freq(i, :), IS_results.cap_np_dt(i, :)', 'Color', Int_colors(i, :), 'Marker', '+', 'MarkerSize', 9, 'LineStyle', '-.', 'LineWidth', 1);
        end
        
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    ax.YScale = 'log'; % for putting the scale in log
    xlim(xlim_array)
    xlabel('Frequency [Hz]');
    ylabel('\omega^{-1} x Im(Z^{-1}) [F/cm^2]');
    legend boxoff
    

