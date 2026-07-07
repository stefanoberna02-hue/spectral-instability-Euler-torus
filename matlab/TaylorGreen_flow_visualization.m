clear; close all; clc;

% Domain
N = 500;
x = linspace(0, 2*pi, N);
y = linspace(0, 2*pi, N);
[X, Y] = meshgrid(x, y);

m = 2;
n = 2;

% Stream function of the (m,n) Taylor--Green vortex
psi = -sin(m*X).*sin(n*Y);

% Velocity field u_E = grad^\perp psi_E = (-psi_y, psi_x)
U = n*sin(m*X).*cos(n*Y);
V = -m*cos(m*X).*sin(n*Y);

% Speed
speed = sqrt(U.^2 + V.^2);

figure('Color','w');
hold on;

% Background: stream function level sets
contourf(X, Y, psi, 35, 'LineColor', 'none');
colormap(parula);
cb = colorbar;
cb.Label.String = '$\psi_E \; (x,y)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;
cb.TickLabelInterpreter = 'latex';

% Streamlines
contour(X, Y, psi, 18, 'k', 'LineWidth', 0.55);

% Vector field, downsampled
step = 0; % set to 0 to remove arrows, other values might be: 100, 50
if step > 0
    quiver(X(1:step:end,1:step:end), ...
           Y(1:step:end,1:step:end), ...
           U(1:step:end,1:step:end), ...
           V(1:step:end,1:step:end), ...
           1.2, 'k', 'LineWidth', 1.0);
end

% Hyperbolic stagnation points:
% sin(mx)=0 and sin(ny)=0
ix = 0:(2*m);
iy = 0:(2*n);
[IX, IY] = meshgrid(ix, iy);

stagX = IX*pi/m;
stagY = IY*pi/n;

plot(stagX(:), stagY(:), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

% Formatting
axis equal tight;
xlim([0 2*pi]);
ylim([0 2*pi]);

xticks([0 pi 2*pi]);
xticklabels({'0','\pi','2\pi'});
yticks([0 pi 2*pi]);
yticklabels({'0','\pi','2\pi'});

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 14);
title('Taylor--Green cellular flow $u_E=\nabla^\perp \; \psi_E$', ...
      'Interpreter', 'latex', 'FontSize', 15);

set(gca, 'FontSize', 12, 'LineWidth', 1);
box on;