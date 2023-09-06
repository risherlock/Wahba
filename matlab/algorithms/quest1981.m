function [C_opt, q_opt] = quest1981(v_b,v_i,w)
%%% Quaternion Estimator (QUEST)
%
% Inputs:
%   v_b = Unit measurement vectors in the spacecraft body frame (3xn).
%   v_i = Corresponding unit vectors known in the inertial frame (3xn).
%   w = Non negative weight assigned to each observations (nx1). 
%
% Outputs:
%   q_opt = Optimal quaternion that transforms vectors in inertial frame to 
%           vectors in body frame ([q0,q1,q2,q3]').
%   C_opt = Optimal rotation matrix that transforms vectors in inertial 
%           frame to vectors in body frame.   
%
% Reference:  
%   [1] Malcolm D. Shuster - Three-axis Attitude Determination from Vector 
%       Observations (1981)
%  
% Rishav (2020-11-6)

% Error tolerance
tolerance = 10e-5;

B = (v_b.*repmat(w,[1,3])')*v_i'; % Eqn.38
Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)]; % Eqn.46
S = B + B'; % Eqn.45
sigma = trace(B); % Eqn.44

% Eqn.63
delta = det(S);
kappa = trace(delta*inv(S));

% Eqn.71
a = sigma^2 - kappa;
b = sigma^2 + Z'*Z;
c = delta + Z'*S*Z;
d = Z'*S^2*Z;
constant = a*b + c*sigma - d;

% Characteristic equation for Newton-Raphson method: (Eqn.70) 
%   f(lambda) = lambda^4 - (a + b)*lambda^2 - c*lambda + constant = 0
%   where, constant = a*b + c*sigma - d
lambda = sum(w); % Shuster1981 and markley1993
last_lambda = 0.0;
while abs(lambda - last_lambda) >= tolerance
    last_lambda = lambda;
    
    f = lambda^4 - (a + b)*lambda^2 - c*lambda + constant;
    f_dot = 4*lambda^3 - 2*(a + b)*lambda - c;
    lambda = lambda - f/f_dot;
end
  
% Eqn.66
omega = lambda;
alpha = omega^2 - sigma^2 + kappa;
beta  = omega - sigma;
gamma = (omega + sigma)*alpha - delta;

% Determine optimal quaternion
X = (alpha*eye(3) + beta*S + S^2)*Z; % Eqn.68
q_opt = [X; gamma]./sqrt(gamma^2 + norm(X)^2); % Eqn.69
C_opt = quaternion_to_dcm(q_opt);
end
