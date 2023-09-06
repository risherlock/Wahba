function [C_opt, q_opt] = yang_manifold2015(v_b,v_i,w) 
%%% Newton's Method on Riemannian manifold
%
% Inputs:
%   v_b = Unit measurement vectors in the spacecraft body frame (3xn).
%   v_i = Corresponding unit vectors known in the inertial frame (3xn).
%
% Output:
%   q_opt = Optimal quaternions that transforms vectors in inertial frame 
%           to vectors in body frame ([q0,q1,q2,q3]').
%   C_opt = Optimal rotation matrix that transforms vectors in inertial 
%           frame to vectors in body frame.
%
% Reference:
%   [1] Yaguang Yang - Attitude Determination using Newton's Method on 
%       Riemannian Manifold (2015)
%
% Rishav (2020-11-8)

% Error tolerance
tolerance = 10e-15;

% Eqn.3
B = (v_b.*repmat(w,[1,3])')*v_i';
Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
K = [B + B'- eye(3)*trace(B), Z; Z', trace(B)];

% Initial values for q_ and lambda.
% Paragraph below Algorithm 2.1  
K_ = -(K - sum(w)*eye(4));
q  = - K_(:,1:3)\K_(:,4);
q_ = [q; 1]/norm([q; 1]); % Initial eigenvector
lambda = q_'*K*q_; % Initial eigenvalue

% K * q_opt = lambda_max * q_opt
% Solving Eqn.6 and Eqn.7 using Newton's method for lambda_max.
while norm(K*q_- lambda*q_) > tolerance % K*q_opt = lambda_max*q_opt
    P_q = eye(4) - q_*q_';
    y   = [(P_q*K*P_q - q_'*K*q_*eye(4)); q_']\[- P_q*K*q_; 0]; % Eqn.6
    q_  = (q_ + y)/norm(q_ + y); % Eqn.7
    lambda = q_'*K*q_; % Eqn.2
end 

% Optimal quaternion
q_opt = q_;
C_opt = quaternion_to_dcm(q_opt);
