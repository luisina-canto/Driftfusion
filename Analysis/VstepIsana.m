function VstepIsana(sol_stepV, figson)

par = sol_stepV.par;
tstep = par.V_fun_arg(3);

tsol = sol_stepV.t;
t = tsol;
Vfun = fun_gen(par.V_fun_type);
Vapp = Vfun(par.V_fun_arg, tsol);
Vappt = Vapp;

V0 = par.V_fun_arg(1);
DeltaV = par.V_fun_arg(2) - V0;

J = dfana.calcJ(sol_stepV);
DeltaJ = J.tot(end, 1) - J.tot(1, 1);
Jt = J.tot(:,1)';

N = length(t);
omega = 2*pi./t;
% Fourier transform
for i = 1:N
    Yreal(i) = (DeltaJ/DeltaV) + (1/DeltaV).*trapz(Jt(1:N-1).*diff(cos(omega(i).*t)));
    Yimag(i) = (1/DeltaV).*trapz(Jt(1:N-1).*diff(sin(omega(i).*t)));
end

Zreal = 1./Yreal;
Zimag = 1./Yimag;

if figson
%     figure(540)
%     semilogx(t, Vappt)
%     xlabel('Time [s]')
%     ylabel('Vapp [V]')
    
    figure(541)
    semilogx(t, J.tot)
    xlabel('Time [s]')
    ylabel('Current density, J [A cm-2]')
    
    figure(542)
    loglog(1./t, (1./(omega.*Zimag)))
    xlabel('Frequency [Hz]')
    ylabel('Zimag-1 omega-1')
    legend('raw', 'smoothed')

    figure(543)
    semilogx(1./t, Zreal)
    xlabel('Frequency [Hz]')
    ylabel('Zreal')
end