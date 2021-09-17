
mubar   = data.mubar;
rho0bar = data.rho0bar;
Hbar    = data.Hbar;
Sbar    = data.Sbar;
CD   = data.CD;
mass = data.mass;

rbar  = state(:, 1);
vbar  = state(:, 2);

hbar = (rbar - 1);
rhobar  = rho0bar * exp(-hbar / Hbar);
qbar    = 0.5*rhobar.*vbar.^2;
Dbar    = qbar.*Sbar.*CD;
gbar = mubar./rbar.^2;


dvdt  = - Dbar./mass - gbar.*sin(fpa);
dvdt_1 = - Dbar./mass;
dvdt_2 = - gbar.*sin(fpa);

figure(11)
hold on
box on
plot(time * ND.TU2s, dvdt.* ((ND.DU2m / ND.TU2s / ND.TU2s) / 1000), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('$\dot{v}$', 'Interpreter', 'LaTeX')
set(gca,'FontSize',16,'FontName','Times');

figure(12)
hold on
box on
plot(time * ND.TU2s, dvdt_1.* ((ND.DU2m / ND.TU2s / ND.TU2s) / 1000), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
% ylabel('$- Dbar./mass$', 'Interpreter', 'LaTeX')
title('$-\frac{D}{m}$',  'Interpreter', 'LaTeX')
set(gca,'FontSize',16,'FontName','Times');

figure(13)
hold on
box on
plot(time * ND.TU2s, dvdt_2.* ((ND.DU2m / ND.TU2s / ND.TU2s) / 1000), 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
% ylabel('$- gbar.*sin(fpa)$', 'Interpreter', 'LaTeX')
title('$-\frac{\mu \sin{\phi}}{r^2}$', 'Interpreter', 'LaTeX')
set(gca,'FontSize',16,'FontName','Times');

figure(14)
hold on
box on
plot(time * ND.TU2s, speed, 'k', 'Linewidth', 1)
xlabel('Time (s)', 'Interpreter', 'LaTeX')
ylabel('Speed (km/s)', 'Interpreter', 'LaTeX')
set(gca,'FontSize',16,'FontName','Times');