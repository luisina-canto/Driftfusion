function Q = getQ(sol)

par = sol.par;
T = par.T;
kB = par.kB;
dev = par.dev_sub;
u = sol.u;
t = sol.t;
x = sol.x;
x_sub = par.x_sub;
pcum1 = par.pcum + 1;   % Cumulative layer points array

%% Get indexes of layers not including interfaces.
% int_index = find(contains(par.layer_type, 'interface'));   % only
% compatible from 2016 onwards
int_logical_layer = strcmp(par.layer_type, 'layer');
int_logical_active = strcmp(par.layer_type, 'active');
int_logical = int_logical_layer + int_logical_active;
loc = find(int_logical); % interface layer locations

%%
rho = dfana.calcrho(sol_OC(4), 'whole');
rho = rho(end,:);

%%
for m = 1:length(loc)
    p_L = pcum1(loc(m)-1);
    p_R = pcum1(loc(m));

end

end