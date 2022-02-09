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

bias = 1e-6;
c = 3e8;

receivers = [250, 250*sqrt(3)/3];

n = 3;

sig_r0 = 10;

sig_particles = 30; % m
sig_time = 1e-9;

velsmoothing = true;
sig_theta = pi/6;
sig_vel = 1;

%% Define Geometry

load traj1

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
    get_noisy_measurements(t, xtrue, pseudolites, 0, sig_r0, 0.00);

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

nparticles = 1000;

% Prior is a uniform 200 x 200 m distribution centered around xtrue(1)
sample_size = [nparticles, 1];

xi0 = xitrue(1) + 200*rand(sample_size) - 100;
eta0 = etatrue(1) + 200*rand(sample_size) - 100;

ss = rel_power(:, 1) * sig_r0;
dists = sqrt(rel_power(:, 1));

b1 = (r(2, 1)*dists(1)-r(1, 1)*dists(2))/(dists(2)-dists(1))/c;
b2 = (r(3, 1)*dists(1)-r(1, 1)*dists(3))/(dists(3)-dists(1))/c;
b3 = (r(3, 1)*dists(2)-r(2, 1)*dists(3))/(dists(3)-dists(2))/c;

b0 = mean([b1, b2, b3]) + ...
    mean(rel_power(:, 1))*std([b1, b2, b3])*randn(sample_size);

particles = [xi0, eta0, b0];

xhat_mmse = zeros(n, T);
xhat_map = zeros(n, T);
Phat = zeros(n, n, T);
% Measurement update and resample step for t = 0

w = getrpdfs(particles, pseudolites, r(:, 1)+b0'*c, ss);
w = w / sum(w);

mu = sum(w.*particles, 1)';
P = wcov(particles, w);

[~, max_ind] = max(w);

xhat_mmse(:, 1) = mu;
xhat_map(:, 1) = particles(max_ind, :);
Phat(:, :, 1) = P;

% Resample particle set
particles = mvnrnd(mu, P, nparticles);

%%

for i = 2:T
    
    if i > 3 && velsmoothing
        % Velocity smoothing dynamics update
        vel = xhat_mmse(:, i-1) - xhat_mmse(:, i-2);
        theta_prev = atan2(vel(2), vel(1));
        thetas = theta_prev + sig_theta*randn(sample_size);
        mags = norm(vel) + sig_vel*randn(sample_size);
        
        particles(:, 1) = particles(:, 1) + mags.*cos(thetas);
        particles(:, 2) = particles(:, 2) + mags.*sin(thetas);
        particles(:, 3) = particles(:, 3) + sig_time*randn(nparticles, 1);
    else
        % random walk dynamics update
        particles(:, 1:2) = particles(:, 1:2) +...
            sig_particles*randn(nparticles, 2);
        particles(:, 3) = particles(:, 3) + sig_time*randn(nparticles, 1);
    end
    
    ss = rel_power(:, i) * sig_r0;
%     disp(ss')

    w = getrpdfs(particles, pseudolites, r(:, i), ss);
    if sum(w) ~= 0
        w = w / sum(w);
    else
        w = ones(sample_size)/nparticles;
    end
    
    mu = sum(w.*particles, 1)';
    P = wcov(particles, w);
    
    [~, max_ind] = max(w);
    
    xhat_mmse(:, i) = mu;
    xhat_map(:, i) = particles(max_ind, :);
    Phat(:, :, i) = P;

    % Resampling from pmf
    particles = sample_pmf(particles, w, nparticles);
    
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
