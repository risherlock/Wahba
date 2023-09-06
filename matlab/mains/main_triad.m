% Test for the implementation of TRIAD method.
% 2020-11-4

clc
clear
close all
format long

% Configuration
bnf = 0.01; % Measurement noise factor

% Generate ground-truth (rotation matrix)
axis = randn(1, 3);
angle = pi * rand();
q_truth = [sin(angle/2) * (axis ./ norm(axis)), cos(angle/2)]';
C_truth = quaternion_to_dcm(q_truth);

% Generate 2 unit inertial vectors
r1 = rand([3, 1]);
r2 = rand([3, 1]);
r1 = r1 ./ norm(r1);
r2 = r2 ./ norm(r2);

% Generate noisy measurement vectors
b1 = C_truth * r1 + bnf * rand(size(r1));
b2 = C_truth * r2 + bnf * rand(size(r2));

% Estimate orientation
[C_hat, q_hat] = triad1964(b1, b2, r1, r2);

% Display result
fprintf('C_truth \n'); disp(C_truth);
fprintf('\nC_hat \n'); disp(C_hat);
fprintf('\nq_hat \n'); disp(q_truth);
fprintf('\nq_hat \n'); disp(q_hat);
