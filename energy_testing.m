%% Input parameters
params_filepath = './Input_files/spiro_mapi_tio2.csv';     % Filepath to the parameters file
output_filename = 'spiro_mapi_tio2';                                    % Filename for output file
light_intensity = 0.5;            % Suns equivalent
Vmax = 1.2;                     % Maximum voltage for cyclic voltammogram
Vmin = -1.2;                    % Minimum voltage for cyclic voltammogram
scan_rate = 1e-3;               % Current-voltage scan rate [Vs-1] 

%% Load in parameters
%par = pc(params_filepath);

%%

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

% EF0_arr = [-4.6, -4.7, -4.8];
% Ncat_arr = [1e15, 1e16, 1e17, 1e18, 1e19];
%EF0_arr = -4.7;
Ncat_arr = [1e16];

% Calculating tghe mobility which varies with the reciprocal of ion
% concentration to give a constant conductivity
mu_carr = 1e19*1e-10./Ncat_arr;

for i = 1:length(Ncat_arr)
    par_struct(i) = par;
    
    par_struct(i).mu_c(3) = mu_carr(i);
    par_struct(i).Ncat = [Ncat_arr(i), Ncat_arr(i), Ncat_arr(i), Ncat_arr(i), Ncat_arr(i)];
    
    
    %par_struct(i).EF0(1) = EF0_arr(i);
    %par_struct(i).Phi_left = EF0_arr(i);
    
    % Everytime you change your parameters in a script use this function:
    par_struct(i) = refresh_device(par_struct(i));
    % soleq = equilibrate(varargin)
    soleq(i) = equilibrate(par_struct(i));
   
end


%% Obtain and plot equilibrium solution
%sol_eq = equilibrate(par, 1);

%% Call function to obtain equilibrium and cyclic voltammogram solutions
sol_CV = doCV(sol_eq(i).el, light_intensity, 0, Vmax, Vmin, scan_rate, 1, 401);

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
figure(2)
dfplot.ELx(sol_CV(2), tplot);
% PLot the carrier densities at time TPLOT
%dfplot.npx(sol_CV, tplot);