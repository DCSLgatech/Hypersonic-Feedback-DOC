function plotSensitivities(time, S, fignum, col)

figure(fignum); 
subplot(1,4,1); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{r H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,2); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{\phi H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,3); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{v H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,4); hold on; grid on;
xlabel('ND Time', 'FontSize', 15);
ylabel('$\bar{S}_{\gamma H_{s}}$','Interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
%
subplot(1,4,1);
plot(time, S(:, 1), 'color', col, 'LineWidth', 1);
%
subplot(1,4,2);
plot(time, S(:, 2), 'color', col, 'LineWIdth', 1);
%
subplot(1,4,3);
plot(time, S(:, 3), 'color', col, 'LineWidth', 1);
%
subplot(1,4,4);
plot(time, S(:, 4), 'color', col, 'LineWIdth', 1);

end