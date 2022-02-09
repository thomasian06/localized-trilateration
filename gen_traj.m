%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ian Thomas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

set(0, 'DefaultAxesLooseInset', [0,0,0,0])
set(0,'defaultAxesFontSize',11)

colors = get(gca, 'colororder');
close all;

%% Declaration of Parameters

% Pseudolite Locations
pseudolites = [0, 0; 5000, 0; 2500, 2500*sqrt(3)];
[npseudolites, dims] =  size(pseudolites);


Dt = 5; % s
% speed = ~1 m/s
xitrue1 = linspace(1250, 2500, 1000);
etatrue1 = linspace(2500*sqrt(3)/2, 0, length(xitrue1));
xitrue2 = linspace(2500, 2500, 1000);
etatrue2 = linspace(0, 2500*sqrt(3)/2, length(xitrue2));
xitrue3 = linspace(2500, 3750, 1000);
etatrue3 = linspace(2500*sqrt(3)/2, 2500*sqrt(3)/2, length(xitrue3));

xitrue = [xitrue1 xitrue2(2:end) xitrue3(2:end)];
etatrue = [etatrue1 etatrue2(2:end) etatrue3(2:end)];

T = length(xitrue);
t = (0:(T-1)) * Dt;
% add process noise of 1 m std
xitrue = xitrue + randn(size(xitrue));
etatrue = etatrue + randn(size(etatrue));

xtrue = [xitrue; etatrue];

figure(1)
hold on;
grid on
axis equal
xlabel('$\xi$ [m]')
ylabel('$\eta$ [m]')
title('Trajectory of Rover')
plot(xitrue, etatrue , 'LineWidth', 1.5);
scatter(pseudolites(1, 1), pseudolites(1, 2), 50, colors(1, :), 'filled')
scatter(pseudolites(2, 1), pseudolites(2, 2), 50, colors(2, :), 'filled')
scatter(pseudolites(3, 1), pseudolites(3, 2), 50, colors(3, :), 'filled')
legend('Path', 'Pseudolite 1', 'Pseudolite 2', 'Pseudolite 3',...
    'location', 'nw')

save traj2 xtrue t T xitrue etatrue Dt pseudolites npseudolites dims
