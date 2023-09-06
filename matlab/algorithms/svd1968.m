function [C_opt, q_opt] = svd1968(v_b,v_i,w)
%%% SVD method
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
%   [1] F. Landis Markley - Attitude Determination using Vector 
%       Observations and Singular Value Decomposition (1968)
%
% Rishav (2020-11-4)

B = (v_b.*repmat(w,[1 3])')*v_i'; % Eqn.5
[U,~,V] = svd(B); % B = U*S*V' (Eqn.6)

d = det(U)*det(V); % Eqn.14
C_opt = U*diag([1, 1, d])*V'; % Eqn.18
q_opt = dcm_to_quaternion(C_opt);
end
