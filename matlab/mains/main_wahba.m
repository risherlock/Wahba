% Test for the implementation of solutions to Wahba's problems.
% 2020-11-4

clc
clear
close all
format long

% Configuration
n = 30; % Number of vector pairs
bnf = 0.01; % Measurement noise factor
w = 0.01 * ones(n, 1); % Measurement weights

% Generate ground-truths (quaternion and DCM)
axis = randn(1, 3);
angle = pi * rand();
q_truth = [sin(angle/2) * (axis./norm(axis)), cos(angle/2)]';
C_truth = quaternion_to_dcm(q_truth);

% Generate n unit inertial vectors, r
r = rand([3, n]);
r = r ./ vecnorm(r);

% Generate noisy measurement vectors, b
b = C_truth * r;
noise = bnf * rand(size(b));
b = b + noise;

% Estimate orientation
[C_hat, q_hat] = quest1981(b, r, w);

% Display result
fprintf('C_truth \n'); disp(C_truth);
fprintf('\nC_hat \n'); disp(C_hat);
fprintf('\nq_hat \n'); disp(q_truth);
fprintf('\nq_hat \n'); disp(q_hat);
