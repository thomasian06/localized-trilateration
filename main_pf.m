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

% Noise Model
R = 20; % m

% Pseudolite Locations
pseudolites = [0, 0; 5000, 0; 2500, 2500*sqrt(3)];
[npseudolites, dims] =  size(pseudolites);

receivers = [250, 250*sqrt(3)/3];

%% Plot Distribution

% Range estimate of 100 meters
nx = 1000;
ny = nx;
xplot = linspace(-300, 300, nx);
yplot = linspace(-300, 300, nx);
[X, Y] = meshgrid(xplot, yplot);

figure(1)
hold on
grid on
% surf(X, Y, rpdf(X, Y, 0, 0, 100, R).*rpdf(X, Y, 150, 0, 100, R));
surf(X, Y, rpdf(X, Y, 0, 0, 250, R));
scatter3(0, 0, 0, 300, 'filled');
xlabel('$\xi^k$ [m]')
ylabel('$\eta^k$ [m]')
zlabel('$p(\xi^k,\eta^k)$')
shading interp
view(-45, 70)

%% Define Geometry

load traj2

figure(2)
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


%% Generate noisy measurements independent of distance from pseudolites

r = reshape(sqrt(sum((xtrue - ...
    reshape(pseudolites', [dims, 1, npseudolites])).^2, 1)),...
    [T, npseudolites])';

% add 20 m std noise to measurements

r = r + 10*randn(size(r));

% [bias, r, rel_power] = ...
%     get_noisy_measurements(t, xtrue, pseudolites, 0, 20, 0.05);

figure(3)
hold on;
grid on;
xlabel('$t$ [s]')
ylabel('$\rho$ [m]')
title('Simulated Range Measurements')
plot(t, r(1, :), 'LineWidth', 1.5)
plot(t, r(2, :), 'LineWidth', 1.5)
plot(t, r(3, :), 'LineWidth', 1.5)
legend('$\rho_1$', '$\rho_2$', '$\rho_3$')


%% Particle Filter

getrpdf = @(particles, pseudolite, r, s) rpdf(particles(:, 1),...
    particles(:, 2), pseudolite(1), pseudolite(2), r, s);
getrpdfs = @(particles, pseudolites, rs,  ss) ...
    getrpdf(particles, pseudolites(1, :), rs(1), ss(1)).*...
    getrpdf(particles, pseudolites(2, :), rs(2), ss(2)).*...
    getrpdf(particles, pseudolites(3, :), rs(3), ss(3));

nparticles = 1000;

% Prior is a uniform 200 x 200 m distribution centered around xtrue(1)
sample_size = [nparticles, 1];

xi0 = xitrue(1) + 200*rand(sample_size) - 100;
eta0 = etatrue(1) + 200*rand(sample_size) - 100;

particles = [xi0, eta0];

ss = [10, 10, 10];

xhat_mmse = zeros(2, T);
xhat_map = zeros(2, T);
Phat = zeros(2, 2, T);
% Measurement update and resample step for t = 0

w = getrpdfs(particles, pseudolites, r(:, 1), ss);
w = w / sum(w);

mu = sum(w.*particles, 1)';
P = wcov(particles, w);

xhat_mmse(:, 1) = mu;
Phat(:, :, 1) = P;

[~, max_ind] = max(w);
xhat_map(:, 1) = particles(max_ind, :);

% Resample particle set
particles = mvnrnd(mu, P, nparticles);


for i = 2:T
    
    % Dynamics update
    particles = particles + 20*randn(nparticles, 2);
    
    % Measurement Update
    w = getrpdfs(particles, pseudolites, r(:, i), ss);
    w = w / sum(w);
    
    [~, max_ind] = max(w);
    xhat_map(:, i) = particles(max_ind, :);
    
    mu = sum(w.*particles, 1)';
    P = wcov(particles, w);
    
    xhat_mmse(:, i) = mu;
    Phat(:, :, i) = P;

    % Resampling (draw from gaussian covariance)
%     fprintf("Timestep %d \n", i)
%     disp(P)
%     fprintf("\n", i)
    particles = mvnrnd(mu, P, nparticles);
    
end

%%

figure(4)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\xi$ [m]')
title('$\xi$ Estimation Error')
plot(t, xhat_mmse(1, :)-xtrue(1, :), 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')

figure(5)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\eta$ [m]')
title('$\eta$ Estimation Error')
plot(t, xhat_mmse(2, :)-xtrue(2, :), 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')


figure(6)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\xi$ [m]')
title('$\xi$ Estimation Error')
plot(t, xhat_map(1, :)-xtrue(1, :), 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')

figure(7)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\eta$ [m]')
title('$\eta$ Estimation Error')
plot(t, xhat_map(2, :)-xtrue(2, :), 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')







%% Functions


function P = wcov(x, w)
    mu = sum(w.*x);
    C = x - mu;
    P = zeros(2, 2);
    for i = 1:length(w)
        P = P + w(i) * C(i, :)' * C(i, :);
    end
end
