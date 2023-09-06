function [C_opt, q_opt] = davenport1968(v_b,v_i,w)
%%% Davenport's q-method
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
%   [1] Paul B. Davenport - A Vector Approach to the Algebra of Rotations 
%       with Applications (1968)
% 
% Rishav (2020-11-4)

% Matrix K (4x4)
B = (v_b.*repmat(w,[1 3])')*v_i';
Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
K = [B + B'- eye(3)*trace(B), Z; Z', trace(B)];

% Eigenvector accociated with largest eigenvalue
[evec, eval] = eig(K);
[max_eval_row,~] = max(eval); 
[~,max_col_index] = max(max_eval_row);

q_opt = evec(:,max_col_index);
C_opt = quaternion_to_dcm(q_opt);
end
