% script for plotting the electrostatic potential of solutions

%will make it general in which variable = the variable we are changing

%% Getting solutions from IS_full_method_script

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
Ncat_arr = [1e14,1e15, 1e16, 1e17, 1e18, 1e19];
<<<<<<< Updated upstream
=======

var=Ncat_arr
>>>>>>> Stashed changes

%% indentify which parameter we want to vary 

% Calculating the mobility which varies with the reciprocal of ion
% concentration to give a constant conductivity
mu_carr = 1e19*1e-10./var;



for i = 1:length(var)
    par_struct(i) = par;
    
    par_struct(i).mu_c(3) = mu_carr(i);
    par_struct(i).Ncat = [var(i), var(i), var(i), var(i), var(i)];
    
    %par_struct(i).EF0(1) = EF0_arr(i);
    %par_struct(i).Phi_left = EF0_arr(i);
    
    % Everytime you change your parameters in a script use this function:
    par_struct(i) = refresh_device(par_struct(i));
    % soleq = equilibrate(varargin)
    soleq(i) = equilibrate(par_struct(i));
end


for i = 1:length(var)
    %% Obtain open circuit initial condition
    % sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
    sol_Rs1e6(i) = lightonRs(soleq(i).ion, 0.5, -100, 1, 1e6, 100);
    sol_OC(i) = RsToClosedCircuit(sol_Rs1e6(i));
end

%% plotting VX

for i = 1:length(var)
                
                df2plot.Vx(length(var),i,var(i),sol_OC(i));
               
end
%% plot recombination rates

for i = 1:length(var)
                
                df2plot.rx(length(var),i,var(i),sol_OC(i));
               
end


%% analyse each zone
%We want to split it up to the 4 boundaries 
% we are interested in 3 voltages: initial, middle, final

%need to find V at x1=200nm and x2=600nm

%find voltage at initial parameter:

V_initial=[]; % create a array to store all the initial voltages
V_final=[] ;% create a array to store all the initial voltages
V_boundary1=[]; %this is the potential between the 1st and 2nd boundary
V_boundary2=[]; %this is the potential between the 3rd and 4th boundary 
inf_pos=[]; %inflextion position of x array corresponding to find V
second_dif=[];
width1=[];
width2=[]


for i = 1:length(Ncat_arr)
    [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_OC(i));
    
    V_ini= V(100,1);
    V_initial(i)=V_ini;
    V_fin=V(100,length(V));
    V_final(i)=V_fin;


    %want the potentials at the boundaries
    %need to find the nearest x values to :
    x1=2e-5;
    x2=6e-5;

    findx1=[];
    findx2=[];

    for j = 1:length(x)
        x_val= x(j);

        %need to find the closest value of x1 fist
        difference1=abs(x_val-x1);
        findx1(j)=difference1;

        %same for x2
        difference2=abs(x_val-x2);
        findx2(j)=difference2;


    end
    %finding the position of the x1 and x2 values
    [M,I]=min(findx1);
    index_x1=I;


    [M,I]=min(findx2);
    index_x2=I;

    %getting the values of V near x1 and x2:

    V_boundary1(i)=V(100,index_x1);
    V_boundary2(i)=V(100,index_x2);


% 
     figure (10+i)
     plot(x(1:740),dif)
%   
% 
%     x_index0=[];
%     x_index01=[];
%     
%     for n=1:length(dif)
%         grad=dif(n);
% 
% %         if grad>0
% %             x_index0(n)=n;
% %             x_position0(n)=x(n);
% %             
% %         elseif 
% %             x(n)
% 
%    
%          if n<=380 %finding the first drop of potential width
%              if grad>0
%                 x_index0(n)=n;
%           
%              end

%          elseif n>380
%              if grad<=-2.5 %second drop of potential
%                 x_index01(n)=n;
%             end                 
%         end
    end    

%now we find the actual position of x
%%
%first width
recom_boundaries=[];
recom_boundaries2=[];

for k=1:length(x_index0)-1
    if abs(x_index0(k+1)-x_index0(k))>1
        recom_boundaries=[recom_boundaries,x(k)];
    end    

end
% 
% width_1=recom_boundaries(2)-recom_boundaries(1);
% width1=[width1, width_1];

for k=380:length(x_index01)-1
    if abs(x_index01(k+1)-x_index01(k))>1
        
        recom_boundaries2=[recom_boundaries2,x(k)];
    end    

end
% 
% width_2=recom_boundaries2(2)-recom_boundaries2(1);
% width2=[width2, width_2];



%%


%% relative recombionation 

%need relative space charge layer widths
%the width is the distance between the 2 curves being flat
%width_r1= width of the 1st region in between green and white 
%width_r2= width of region between white and blue 



%%

%recombination rates

%% integrate rx to give tot_r (area under rx graph)

x=[];
y=[];
tot_r=[];

for i=1:length(var)
    [x y]=rx(length(var),i,var(i),sol_OC(i));
   tot_r(i)=trapz(x*1e-7,y);

end   

%% finding recombination for chosen regions

%setting up the values 

Ncat=[0 0 0 0 0 0; 1.68E-05	2.00E-05	4.77E-05	4.77E-05	6.00E-05	6.80E-05
;1.55E-05	2.00E-05	4.50E-05	4.50E-05	6.00E-05	6.36E-05
;1.04969E-05	2.00E-05	0.000039768	4.17544E-05	6.00E-05	6.85828E-05
;0.000010249	2.00E-05	0.000229958	5.75298E-05	6.00E-05	6.39767E-05;...
1.02E-05	2.00E-05	2.52E-05	0.0000581	6.00E-05	0.0000647];

%Ncat=[{},{1.68E-05,2.00E-05,4.77E-05,4.77E-05,6.00E-05,6.80E-05},{1.55E-05,2.00E-05,4.50E-05,4.50E-05,6.00E-05,6.36E-05},{1.04969E-05,2.00E-05,0.000039768,4.17544E-05,6.00E-05,6.85828E-05}...
%,{0.000010249,2.00E-05,0.000229958,5.75298E-05,6.00E-05,6.39767E-05},{1.02E-05,2.00E-05,2.52E-05,0.0000581,6.00E-05,0.0000647}];

%%



for i=1:length(var)
   [x y]=rx(length(var),i,var(i),sol_OC(i));
   x1=x(Ncat(i,1):Ncat(i,2));
   y1=y(Ncat(i,1):Ncat(i,2));


end   



 

%%



%%

 function [x y]=rx(leng,loop_index,value,varargin)
            % Recombination rates as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            x_sub = par.x_sub;
            r = dfana.calcr(sol, "sub");
            y=r.tot(100,:);
            x=x_sub*1e7;
            
           
        end

%% funciton to get the inflection point 
%for electrostatic graphs  which dont have a clear straight line at the mid
%need to make a function to find 2nd inflection point

function curv=inflection(x,V)
    Vp=V(100,:);
    dVp=Vp(2:end)-Vp(1:end-1);
    dx=x(2:end)-x(1:end-1);
    d2Vp=dVp(2:end)-dVp(1:end-1);
    d2x=dx(2:end)-dx(1:end-1);
    curv=d2Vp./d2x;
    
end



%% 

function dif=differential(x,V)
    Vp=V(100,:);
    dVp=Vp(2:end)-Vp(1:end-1);
    dx=x(2:end)-x(1:end-1);
    dif=dVp./dx;

end





