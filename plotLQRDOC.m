function plotLQRDOC(lqrdoc_solution, ND, C)

t = lqrdoc_solution.phase.time * ND.TU2s;
state = lqrdoc_solution.phase.state;
control = lqrdoc_solution.phase.control;

alt       = (state(:,1) * ND.DU2m - C.Rm)/1000;
speed     = state(:,2)* (ND.DU2m / ND.TU2s) / 1000;
fpa       = state(:,3)*180/pi;
dran      = state(:,4) * ND.DU2m / 1000;
xplot     = [alt, speed, fpa, dran];
SHplot    = abs([state(:,5), state(:,6)/ND.TU2s, state(:,7)/ND.DU2m, state(:,8)]); 


    % plotting
    figure(1)
    subplot(2,2,1)
    hold on
    plot(t, xplot(:,1), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Altitude (km)', 'Interpreter', 'LaTeX')
    
    subplot(2,2,2)
    hold on
    plot(t, xplot(:,2), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Speed (km/s)', 'Interpreter', 'LaTeX')
    
    subplot(2,2,3)
    hold on
    plot(t, xplot(:,3), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Flight path angle (deg)', 'Interpreter', 'LaTeX')
    
    subplot(2,2,4)
    hold on
    plot(t, xplot(:,4), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('Downrange (km)', 'Interpreter', 'LaTeX')
    
    figure(2)
    hold on
    plot(t, control, '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$\cos \sigma$', 'Interpreter', 'LaTeX')
    legend('Open-loop', 'LQR-DOC Coupled')
    
    figure(3)
    subplot(2,2,1)
    hold on
    plot(t, SHplot(:,1), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_1$', 'Interpreter', 'LaTeX')
    legend('Open-loop', 'LQR-DOC Decoupled', 'LQR-DOC Coupled')
    
    subplot(2,2,2)
    hold on
    plot(t, SHplot(:,2), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_2$', 'Interpreter', 'LaTeX')
    
    subplot(2,2,3)
    hold on
    plot(t, SHplot(:,3), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_3$', 'Interpreter', 'LaTeX')
    
    subplot(2,2,4)
    hold on
    plot(t, SHplot(:,4), '--', 'color', 'r', 'Linewidth', 1)
    xlabel('Time (s)', 'Interpreter', 'LaTeX')
    ylabel('$S_4$', 'Interpreter', 'LaTeX')

end