function [ph, x, y] = conf_ellipse(m, P, N, conf, plot_flag, varargin)

t = linspace(0, 2 * pi, N);
R = chol(inv(P));
z = sqrt(-2 * log(1 - conf)) * [cos(t); sin(t)];
X = (R \ z);

x = m(1) + X(1, :);
y = m(2) + X(2, :);

if plot_flag
    ph = plot(x, y, varargin{:});
end
ph = 0;

end