function plotSensitivities(ND, time, S, fignum, col)

figure(fignum); 
subplot(4,1,1); hold on; grid on;
xlabel('Time (min)', 'FontSize', 15);
ylabel('$\bar{S}_{r H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,1,2); hold on; grid on;
xlabel('Time (min)', 'FontSize', 15);
ylabel('$\bar{S}_{\phi H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,1,3); hold on; grid on;
xlabel('Time (min)', 'FontSize', 15);
ylabel('$\bar{S}_{v H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,1,4); hold on; grid on;
xlabel('Time (min)', 'FontSize', 15);
ylabel('$\bar{S}_{\gamma H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(4,1,1);
plot(time*ND.TU2s/60, S(:, 1), 'color', col, 'LineWidth', 1);
%
subplot(4,1,2);
plot(time*ND.TU2s/60, S(:, 2), 'color', col, 'LineWIdth', 1);
%
subplot(4,1,3);
plot(time*ND.TU2s/60, S(:, 3), 'color', col, 'LineWidth', 1);
%
subplot(4,1,4);
plot(time*ND.TU2s/60, S(:, 4), 'color', col, 'LineWIdth', 1);

end