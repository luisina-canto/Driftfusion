function [Zreal, Zimag, Zmag, Zphase] = VstepISana(sol_stepV, figson)
% F_CUTOFF = Cutoff frequency for plots
tic
disp('Starting Vstep IS analysis')
par = sol_stepV.par;

tsol = sol_stepV.t;
t = tsol;               
t = t(3:end) - t(3);                    % Remove initial ramp
f_cutoff = 1/(10^(log10(tsol(2)) + 2));      % High frequency cutoff for plots

V0 = par.V_fun_arg(1);
DeltaV = par.V_fun_arg(2) - V0;

J = dfana.calcJ(sol_stepV);
Jt = (J.tot(:, 1))';
Jt = Jt(3:end);                 % Remove initial ramp
% zeroing
Jt = abs(Jt - Jt(end));
DeltaJ = abs(Jt(end) - Jt(1));

N = length(t);
omega = 2*pi./t;
% Fourier transform method from Neukom 2019 Thesis, Section 8.1.1
for i = 1:N
    Yreal(i) = (DeltaJ/DeltaV) + (1/DeltaV).*sum(Jt(1:N-1).*(abs(cos(omega(i).*t(1:N-1))) - abs(cos(omega(i).*t(2:N)))));
    Yimag(i) = (1/DeltaV).*sum(Jt(1:N-1).*(abs(sin(omega(i).*t(2:N))) - abs(sin(omega(i).*t(1:N-1)))));
end

Zreal = 1./Yreal;
Zimag = 1./Yimag;
Zmag = (Zreal.^2 + Zimag.^2).^0.5;
Zphase = (rad2deg(atan2(Zimag,Zreal)));

if figson
    figure(51)
    loglog(t, Jt)
    xlabel('Time [s]')
    ylabel('Current density, J [A cm-2]')
    xlim([t(1), t(end)])
    
    figure(52)
    loglog(1./t, (1./abs((omega.*-Zimag))))
    xlabel('Frequency [Hz]')
    ylabel('imag(Z^{-1}) \omega^{-1} (F cm-2)')
    xlim([1/t(end), f_cutoff])

    figure(53)
    loglog(1./t, Zmag, '-')
    xlabel('Frequency [Hz]')
    ylabel('|Z| (\Omega cm2)')
    xlim([1/t(end), f_cutoff])
    
    figure(54)
    semilogx(1./t, Zphase)
    xlabel('Frequency [Hz]')
    ylabel('Phase (degrees)')
    xlim([1/t(end), f_cutoff])
    %ylim([-180, 10])
    
    figure(55)
    plot(Zreal(t > 1./f_cutoff), Zimag(t > 1./f_cutoff))
    xlabel('Zreal (\Omega cm2)')
    ylabel('Zimag (\Omega cm2)')
end

disp('Vstep IS analysis complete')
toc
end