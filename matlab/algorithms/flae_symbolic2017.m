function [C_opt, q_opt] = flae_symbolic2017(v_b, v_i, w)
%%% Fast Linear Quaternion Attitude Estimator (FLAE) using symbolic method
%
% Inputs:
%   v_b = Unit measurement vectors in the spacecraft body frame (3xn).
%   v_i = Corresponding unit vectors known in the inertial frame (3xn).
%   w = Non negative weight assigned to each observations (nx1). 
%
% Outputs:
%   C_opt = Optimal rotation matrix that transforms vectors in inertial 
%           frame to vectors in body frame.   
%   q_opt = Optimal quaternion that transforms vectors in inertial frame to
%           vectors in body frame ([q1,q2,q3,q0]').
%
% Reference:  
%   [1] Jin Wu (et al.) - Fast Linear Quaternion Attitude Estimator Using 
%       Vector Observations (2017)
% 
% Rishav (2020-12-20)

% Hx, Hy, Hz (1x3 matrices) (Eqn.42)
B = (v_b.*repmat(w,[1,3])')*v_i';
H_x = B(:,1)';
H_y = B(:,2)';
H_z = B(:,3)';

% Eqn.45
W = zeros(4);
W(1,1) = H_x(1) + H_y(2) + H_z(3);
W(1,2) = - H_y(3) + H_z(2);
W(1,3) = - H_z(1) + H_x(3);
W(1,4) = - H_x(2) + H_y(1);
W(2,1) = - H_y(3) + H_z(2);
W(2,2) = H_x(1) - H_y(2) - H_z(3);
W(2,3) = H_x(2) + H_y(1);
W(2,4) = H_x(3) + H_z(1);
W(3,1) = - H_z(1) + H_x(3);
W(3,2) = H_x(2) + H_y(1);
W(3,3) = H_y(2) - H_x(1) - H_z(3);
W(3,4) = H_y(3) + H_z(2);
W(4,1) = -H_x(2) + H_y(1);
W(4,2) = H_x(3) + H_z(1);
W(4,3) = H_y(3) + H_z(2);
W(4,4) = H_z(3) - H_y(2) - H_x(1);

% Eqn.49
tau_1 = -2 * ( H_x(1)^2 + H_x(2)^2 + H_x(3)^2 ...
    + H_y(1)^2 + H_y(2)^2 + H_y(3)^2 ...
    + H_z(1)^2 + H_z(2)^2 + H_z(3)^2 );
tau_2 = 8 * ( H_x(3) * H_y(2) * H_z(1) - H_x(2) * H_y(3) * H_z(1) - ...
    H_x(3) * H_y(1) * H_z(2) + H_x(1) * H_y(3) * H_z(2) + ...
    H_x(2) * H_y(1) * H_z(3) - H_x(1) * H_y(2) * H_z(3) );
tau_3 = det(W);

% Symbolic approach
T0 = 2*tau_1^3 + 27*tau_2^2 - 72*tau_1*tau_3;
T1 = (T0 + sqrt( - 4*(tau_1^2 + 12*tau_3)^3 + T0^2))^(1/3);
T2 = sqrt( - 4 * tau_1 + 2^(4/3)*(tau_1^2 + 12*tau_3)/T1 + 2^(2/3)*T1);

% Eqn.53
sqrt6x12    = 29.393876913398135;
sqrt6x2_inv = 0.20412414523193150818310700622549;
lambda1 =   sqrt6x2_inv*(T2 - sqrt(-T2^2 - 12*tau_1 - 12*sqrt6x12*tau_2/T2));
lambda2 =   sqrt6x2_inv*(T2 + sqrt(-T2^2 - 12*tau_1 - 12*sqrt6x12*tau_2/T2));
lambda3 = - sqrt6x2_inv*(T2 + sqrt(-T2^2 - 12*tau_1 + 12*sqrt6x12*tau_2/T2));
lambda4 = - sqrt6x2_inv*(T2 - sqrt(-T2^2 - 12*tau_1 + 12*sqrt6x12*tau_2/T2));
lambda = max(real([lambda1, lambda2, lambda3, lambda4]));
N = W - lambda * eye(4); % Eqn.54

% Transform N to N' via elementary row operations (Eqn.55)
% N -> N' = [ 1,  0,  0,  chi;
%             0,  1,  0,  rho;
%             0,  0,  1,  nu;
%             0,  0,  0,  zeta ]
% zeta is around 1x10-15 (i.e. approx zero)
N(1,:) = N(1,:) / N(1,1);
N(2,:) = N(2,:) - N(2,1)*N(1,:);
N(3,:) = N(3,:) - N(3,1)*N(1,:);
N(4,:) = N(4,:) - N(4,1)*N(1,:);

N(2,:) = N(2,:) / N(2, 2);
N(1,:) = N(1,:) - N(1,2)*N(2,:);
N(3,:) = N(3,:) - N(3,2)*N(2,:);
N(4,:) = N(4,:) - N(4,2)*N(2,:);

N(3,:) = N(3,:) / N(3, 3);
N(1,:) = N(1,:) - N(1,3)*N(3,:);
N(2,:) = N(2,:) - N(2,3)*N(3,:);
N(4,:) = N(4,:) - N(4,3)*N(3,:);

% Fundamental solutions: chi, rho and nu (Eqn.57) 
chi = N(1,4);
rho = N(2,4);
nu  = N(3,4);
q = [rho, nu, -1, chi]'; % Eqn.58

% Optimal unit quaternion
q_opt = q / norm(q);
C_opt = quaternion_to_dcm(q_opt);
end
