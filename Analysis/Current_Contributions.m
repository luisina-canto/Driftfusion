function Jcomp = current_contributions(sol, figson)
% Script to plot the currents from radiative and non-radiative (split into SRH
% and VSR) losses, alongside the current measured in a JV sweep and the
% generation current.
% JCOMP is a structure containing the different current outputs as a
% function of time:
% J_GEN = Generation current
% J_BTB = Band-to-band (radiative) recombination current
% J_SRH_BULK = Shockley-Read-Hall bulk recombiantion current
% J_VSR = Shockley-Read-Hal internal interfaces recombination currents. The
% row number defines the interface numbered from left to right
% J_SURF_L = Left-hand boundary SRH surface recombination current
% J_SURF_R = Right-hand boundary SRH surface recombination current
% J_TOT = Total current

par = sol.par;
t = sol.t;
pcum1 = par.pcum + 1;   % Cumulative layer points array
%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr and J_ext
num_values = length(sol.t);
e = sol.par.e;
loss_currents = dfana.calcr(sol, 'sub');
x = sol.par.x_sub;
gxt = dfana.calcg(sol);
J = dfana.calcJ(sol);
j_surf_rec = dfana.calcj_surf_rec(sol);
V = dfana.calcVapp(sol);

%% Get indexes of interfaces
%int_index = find(contains(par.layer_type, 'interface'));   % only
%compatible from 2016 onwards
int_logical = zeros(1, length(par.layer_type));
for i = 1:length(par.layer_type)
    int_logical(i) = any(any(strcmp(par.layer_type(i), {'interface', 'junction'})));
end
loc = find(int_logical); % interface layer locations

J_vsr = zeros(length(loc), length(t));
for i = 1:length(loc)   
    % Check location of interfacial surface carrier densities
    p_L = pcum1(loc(i)-1);
    p_R = pcum1(loc(i));

    J_vsr(i,:) = e*trapz(x(p_L:p_R), loss_currents.vsr(:, p_L:p_R), 2)';
end
%forward sweep
J_gen = -e*trapz(x, gxt, 2)';
J_btb = e*trapz(x, loss_currents.btb, 2)';
J_srh_bulk = e*trapz(x, loss_currents.srh, 2)';
J_surf_l = e*(j_surf_rec.l)';
J_surf_r = e*(j_surf_rec.r)';
J_tot = J.tot(:,1)';

% Package into output structure
Jcomp.J_gen = J_gen;
Jcomp.J_btb = J_btb;
Jcomp.J_srh_bulk = J_srh_bulk;
Jcomp.J_surf_l = J_surf_l;
Jcomp.J_surf_r = J_surf_r;
Jcomp.J_vsr = J_vsr;
Jcomp .J_tot = J_tot;

if figson
    %% Plot contributons to the current
    figure(300)
    plot(V, J_gen, '--', V, J_btb, V, J_srh_bulk, V, J_surf_l, V, J_surf_r)
    hold on
    for i = 1:length(loc)
        plot(V, J_vsr(i, :))
    end
    plot(V, J_tot, 'k')
    
    plot(V(1:num_values), zeros(1,num_values), 'black', 'LineWidth', 1)
    hold off
    xlabel('Voltage (V)')
    ylabel('Current Density (Acm^{-2})')
    legend({'J_{gen}', 'J_{rad}', 'J_{SRH, bulk}', 'J_{surf, l}', 'J_{surf, r}', 'J_{interface 1}', 'J_{interface 2}', 'J_{ext}'}, 'Location', 'bestoutside')

    figure(301)
    plot(t, J_gen, '--', t, J_btb, t, J_srh_bulk, t, J_surf_l, t, J_surf_r)
    hold on
    for i = 1:length(loc)
        plot(t, J_vsr(i, :))
    end
    plot(t, J_tot, 'k')
    
    plot(t, zeros(1, num_values), 'black', 'LineWidth', 1)
    hold off
    xlim([0, t(end)])
    xlabel('Times (s)')
    ylabel('Current Density (Acm^{-2})')
    legend({'J_{gen}', 'J_{rad}', 'J_{SRH, bulk}', 'J_{surf, l}', 'J_{surf, r}', 'J_{interface 1}', 'J_{interface 2}', 'J_{ext}'}, 'Location', 'bestoutside')

end

end