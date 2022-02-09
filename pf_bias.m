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
set(0,'defaultAxesFontSize',12)

colors = get(gca, 'colororder');
close all;

rng(101)

%% Declaration of Parameters

% Pseudolite Locations
pseudolites = [0, 0; 5000, 0; 2500, 2500*sqrt(3)];
[npseudolites, dims] =  size(pseudolites);

cpoly = [0];
% cpoly = [3.5556e-14, 5.3333e-10, 1e-6];
c = 3e8;

n = 3;

sig_r0 = 10;
dfac = 0.05;

sig_particles = 30; % m
sig_time = 1e-7;

velsmoothing = false;
sig_theta = pi/6;
sig_vel = 20;

nls_cost = @(x, p, r, s) sum((1./s.^2).*...
    (sqrt((x(1)-p(:, 1)).^2+(x(2)-p(:, 2)).^2)-r+c*x(3)).^2);


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

figure(20)
hold on
grid on
axis equal
xlabel('$t$ [s]')
ylabel('$x$ [m]')
title('State vs $t$')
plot(t, xtrue , 'LineWidth', 1.5);
legend('$\xi$', '$\eta$', 'location', 'se')


%% Generate noisy measurements independent of distance from pseudolites

[bias, r, rel_sig] = ...
    get_noisy_measurements(t, xtrue, pseudolites, cpoly, sig_r0, dfac);

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

ss = rel_sig(:, 1) * sig_r0;
dists = rel_sig(:, 1);

b1 = (r(2, 1)*dists(1)-r(1, 1)*dists(2))/(dists(2)-dists(1))/c;
b2 = (r(3, 1)*dists(1)-r(1, 1)*dists(3))/(dists(3)-dists(1))/c;
b3 = (r(3, 1)*dists(2)-r(2, 1)*dists(3))/(dists(3)-dists(2))/c;

b0 = mean([b1, b2, b3]) + ...
    mean(rel_sig(:, 1))*std([b1, b2, b3])*randn(sample_size);

particles = [xi0, eta0, b0];

xhat_NLS = zeros(n, T);
xhat_mmse = zeros(n, T);
xhat_map = zeros(n, T);
Phat = zeros(n, n, T);
% Measurement update and resample step for t = 0

% w = getrpdfs(particles, pseudolites, r(:, 1)+b0'*c, ss);
w = getrpdfs(particles, pseudolites, r(:, 1), ss);
if sum(w) ~= 0
    w = w / sum(w);
else
    w = ones(sample_size)/nparticles;
end

Ness = zeros(1, T);
Ness(1) = 1/sum(w.^2);

mu = sum(w.*particles, 1)';
P = wcov(particles, w);

% xhat_NLS(:, 1) = 
xhat_mmse(:, 1) = mu;
Phat(:, :, 1) = P;

[~, max_ind] = max(w);
xhat_map(:, 1) = particles(max_ind, :);

xhat_NLS(:, 1) = fminsearch(@(x) nls_cost(x, pseudolites, r(:, 1), ss), [1250;0;0]);

% Resample particle set
particles = sample_pmf(particles, w, nparticles);

%%

for i = 2:T

    disp(i)
    if i > 3 && velsmoothing
        % Velocity smoothing dynamics update
        vel = xhat_mmse(1:2, i-1) - xhat_mmse(1:2, i-2);
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
    
    ss = (rel_sig(:, i)) * sig_r0;
    dists = rel_sig(:, i);
    
    % Estimate clock bias roughly from relative power
    b1 = (r(2, i)*dists(1)-r(1, i)*dists(2))/(dists(2)-dists(1))/c;
    b2 = (r(3, i)*dists(1)-r(1, i)*dists(3))/(dists(3)-dists(1))/c;
    b3 = (r(3, i)*dists(2)-r(2, i)*dists(3))/(dists(3)-dists(2))/c;
    
    mb = mean([b1 b2 b3]);
    stdb = std([b1 b2 b3]);
    particles(:, 3) = mb + stdb*randn(sample_size);
    % Measurement Update
%     w = getrpdfs(particles, pseudolites, r(:, i)+particles(:, 3)'*c, ss);
    w = getrpdfs(particles, pseudolites, r(:, i), ss);
    if sum(w) ~= 0
        w = w / sum(w);
    else
        w = ones(sample_size)/nparticles;
    end
    Ness(i) = 1/sum(w.^2);
%     w = w .* normpdf(particles(:, 3), mb,...
%         mean(rel_sig(:, i))*stdb);
%     if sum(w) ~= 0
%         w = w / sum(w);
%     else
%         w = ones(sample_size)/nparticles;
%     end
    
    [~, max_ind] = max(w);
    xhat_map(:, i) = particles(max_ind, :);
    
    mu = sum(w.*particles, 1)';
    P = wcov(particles, w);
    
    xhat_mmse(:, i) = mu;
    Phat(:, :, i) = P;

    xhat_NLS(:, i) = fminsearch(@(x)...
        nls_cost(x, pseudolites, r(:, i), ss), xhat_NLS(:, i-1));
    
    % Resampling (draw from gaussian covariance)
    particles = sample_pmf(particles, w, nparticles);
    
end

%%

figure(4)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\xi$ [m]')
title('$\xi$ MMSE Estimation Error')
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
title('$\eta$ MMSE Estimation Error')
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
ylabel('$e_b$ [s]')
title('$b$ MMSE Estimation Error')
plot(t, xhat_mmse(3, :)-bias, 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(3, 3, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(3, 3, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')


figure(7)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\xi$ [m]')
title('$\xi$ MAP Estimation Error')
plot(t, xhat_map(1, :)-xtrue(1, :), 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')

figure(8)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\eta$ [m]')
title('$\eta$ MAP Estimation Error')
plot(t, xhat_map(2, :)-xtrue(2, :), 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')

figure(9)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_b$ [s]')
title('$b$ MAP Estimation Error')
plot(t, xhat_map(3, :)-bias, 'LineWidth', 1.5);
plot(t, 2*sqrt(reshape(Phat(3, 3, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
plot(t, -2*sqrt(reshape(Phat(3, 3, :), 1, [])), '--', 'LineWidth', 1.5,...
    'Color', colors(2, :), 'HandleVisibility', 'off')
legend('Error', '2$\sigma$', 'location', 'se')

figure(10)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\xi$ [m]')
title('$\xi$ NLS Estimation Error')
plot(t, xhat_NLS(1, :)-xtrue(1, :), 'LineWidth', 1.5);
% plot(t, 2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
%     'Color', colors(2, :))
% plot(t, -2*sqrt(reshape(Phat(1, 1, :), 1, [])), '--', 'LineWidth', 1.5,...
%     'Color', colors(2, :), 'HandleVisibility', 'off')
% legend('Error', '2$\sigma$', 'location', 'se')

figure(11)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_\eta$ [m]')
title('$\eta$ NLS Estimation Error')
plot(t, xhat_NLS(2, :)-xtrue(2, :), 'LineWidth', 1.5);
% plot(t, 2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
%     'Color', colors(2, :))
% plot(t, -2*sqrt(reshape(Phat(2, 2, :), 1, [])), '--', 'LineWidth', 1.5,...
%     'Color', colors(2, :), 'HandleVisibility', 'off')
% legend('Error', '2$\sigma$', 'location', 'se')

figure(12)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$e_b$ [s]')
title('$b$ NLS Estimation Error')
plot(t, xhat_NLS(3, :)-bias, 'LineWidth', 1.5);
% plot(t, 2*sqrt(reshape(Phat(3, 3, :), 1, [])), '--', 'LineWidth', 1.5,...
%     'Color', colors(2, :))
% plot(t, -2*sqrt(reshape(Phat(3, 3, :), 1, [])), '--', 'LineWidth', 1.5,...
%     'Color', colors(2, :), 'HandleVisibility', 'off')
% legend('Error', '2$\sigma$', 'location', 'se')


figure(13)
hold on
grid on
xlabel('$t$ [s]')
ylabel('$N_\mathrm{ess}$')
title('Effective Sample Size vs. time')
plot(t, Ness, 'LineWidth', 1.5);


%% Functions


function P = wcov(x, w)
    mu = sum(w.*x);
    C = x - mu;
    P = zeros(3, 3);
    for i = 1:length(w)
        P = P + w(i) * C(i, :)' * C(i, :);
    end
end
