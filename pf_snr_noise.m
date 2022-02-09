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

bias = 1e-6;
c = 3e8;

receivers = [250, 250*sqrt(3)/3];

n = 3;

sig_r0 = 10;

%% Plot Distribution

% Range estimate of 100 meters
nx = 100;
ny = nx;
xplot = linspace(-150, 150, nx);
yplot = linspace(-150, 150, nx);
[X, Y] = meshgrid(xplot, yplot);

% rpdf = @(X, Y, x, y, r, s) exp(-0.5*((sqrt((X-x).^2+(Y-y).^2)-r)./s).^2);

figure(1)
hold on
grid on
surf(X, Y, rpdf(X, Y, 0, 0, 100, R).*rpdf(X, Y, 150, 0, 100, R));
xlabel('$x$ [m]')
ylabel('$y$ [m]')
zlabel('$p(x,y)$')
shading interp
view(-45, 58)

%% Define Geometry

% start at -5000, 0 then walk to pseudolite 1, then walk 45 degs north
Dt = 5; % s
% speed = ~1 m/s
xitrue1 = -5000:5:0;
etatrue1 = linspace(0, 50, length(xitrue1));
xitrue2 = (0:5:10000) / sqrt(2);
etatrue2 = linspace(0, 10000, length(xitrue2)) / sqrt(2) + 50;

xitrue = [xitrue1 xitrue2(2:end)];
etatrue = [etatrue1 etatrue2(2:end)];
biastrue = bias*ones(size(xitrue));

T = length(xitrue);
t = (0:(T-1)) * Dt;
% add process noise of 1 m std
xitrue = xitrue + randn(size(xitrue));
etatrue = etatrue + randn(size(etatrue));

xtrue = [xitrue; etatrue];

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

[bias, r, rel_power] = ...
    get_noisy_measurements(t, xtrue, pseudolites, 0, sig_r0, 0.05);

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
    getrpdf(particles, pseudolites(1, :), rs(1, :), ss(1)).*...
    getrpdf(particles, pseudolites(2, :), rs(2, :), ss(2)).*...
    getrpdf(particles, pseudolites(3, :), rs(3, :), ss(3));

nparticles = 100;

% Prior is a uniform 200 x 200 m distribution centered around xtrue(1)
sample_size = [nparticles, 1];

xi0 = xitrue(1) + 200*rand(sample_size) - 100;
eta0 = etatrue(1) + 200*rand(sample_size) - 100;
b0 = bias(1)*ones(sample_size);

particles = [xi0, eta0, b0];

ss = rel_power(:, 1) * sig_r0;

xhat_mmse = zeros(n, T);
xhat_map = zeros(n, T);
Phat = zeros(n, n, T);
% Measurement update and resample step for t = 0

w = getrpdfs(particles, pseudolites, r(:, 1)-b0(1)*c, ss);
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
    
    disp(i)
    % Dynamics update
    particles(:, 1:2) = particles(:, 1:2) + 10*randn(nparticles, 2);
    particles(:, 3) = particles(:, 3) + 0*randn(nparticles, 1);
    
    ss = rel_power(:, i) * sig_r0;
    dists = sqrt(1./rel_power);
    
    % Estimate clock bias roughly from relative power
%     c = 
    
    % Measurement Update
    w = getrpdfs(particles, pseudolites, r(:, i)-particles(:, 3)'*c, ss);
    w = w / sum(w);
    
    mu = sum(w.*particles, 1)';
    P = wcov(particles, w);
    
    [~, max_ind] = max(w);
    xhat_map(:, i) = particles(max_ind, :);
    
    xhat_mmse(:, i) = mu;
    Phat(:, :, i) = P;

    % Resampling (draw from gaussian covariance)
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
    P = zeros(3, 3);
    for i = 1:length(w)
        P = P + w(i) * C(i, :)' * C(i, :);
    end
end
