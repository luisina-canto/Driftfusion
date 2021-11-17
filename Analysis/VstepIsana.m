function VstepIsana(sol_stepV, figson)

par = sol_stepV.par;
tstep = par.V_fun_arg(3);

tsol = sol_stepV.t;
t = tsol(tsol > tstep); % Remove initial ramp
t = t - t(1);           % Zeroing
Vfun = fun_gen(par.V_fun_type);
Vapp = Vfun(par.V_fun_arg, tsol);
Vappt = Vapp(tsol > tstep);

V0 = par.V_fun_arg(1);
DeltaV = par.V_fun_arg(2) - V0;

J = dfana.calcJ(sol_stepV);
DeltaJ = abs(J.tot(end, 1) - J.tot(1, 1));
Jt = (J.tot((tsol > tstep), 1))';
% zeroing
Jt = (Jt - Jt(end));

N = length(t);
omega = 2*pi./t;
% Fourier transform - broken into time units 
for i = 1:N
    Yreal(i) = (DeltaJ/DeltaV) + (1/DeltaV).*sum(Jt(1:N-1).*(cos(omega(i).*t(1:N-1)) - cos(omega(i).*t(2:N))));
    Yimag(i) = (1/DeltaV).*sum(Jt(1:N-1).*(sin(omega(i).*t(2:N)) - sin(omega(i).*t(1:N-1))));
end

Zreal = 1./Yreal;
Zimag = 1./Yimag;

if figson
%     figure(50)
%     semilogx(t, Vappt)
%     xlabel('Time [s]')
%     ylabel('Vapp [V]')
    
    figure(51)
    loglog(t, Jt)
    xlabel('Time [s]')
    ylabel('Current density, J [A cm-2]')
    
    figure(52)
    loglog(1./t, (1./(omega.*Zimag)))
    xlabel('Frequency [Hz]')
    ylabel('imag(Z^{-1}) \omega^{-1}')

    figure(53)
    semilogx(1./t, Zreal)
    xlabel('Frequency [Hz]')
    ylabel('real(Z)')
end