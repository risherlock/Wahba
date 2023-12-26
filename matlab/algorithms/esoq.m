function [C, q] = esoq(vb, vi, w)
%%% Estimator of the Optimal Quaternion (ESOQ)
%
% Inputs:
%   vb = Unit measurement vectors in the spacecraft body frame (3xn).
%   vi = Corresponding unit vectors known in the inertial frame (3xn).
%   w  = Non negative weight assigned to each observations (nx1).
%
% Outputs:
%   q = Optimal quaternion that transforms vectors in inertial frame to
%           vectors in body frame ([q0,q1,q2,q3]').
%   C = Optimal rotation matrix that transforms vectors in inertial
%           frame to vectors in body frame.
%
% Reference:
%   [1] Mortari - ESOQ, A Closed-Form Solution to the Wahba Problem (1997)
%
% Rishav (2023-12-24)

  % Basic matrices and vectors
  B = (vb .* repmat(w, [1,3])') * vi'; % Eqn.4
  z = [B(2,3) - B(3,2); B(3,1) - B(1,3); B(1,2) - B(2,1)];
  S = B + B'- eye(3) * trace(B); % Eqn.6
  K = [S, z; z', trace(B)];

  % Largest eigenvalue of K, lambda
  adj_BBt = det(B + B') * inv(B + B');
  adj_K = det(K) * inv(K);
  b = -2 * trace(B)^2 + trace(adj_BBt) - z' * z; % Eqn.8
  c = -trace(adj_K);
  d = det(K);
  p = (b/3)^2 + 4*d/3; % Eqn.11
  q_ = (b/3)^3 - 4*d*b/3 + c^2/2;
  u1 = 2 * sqrt(p) * cos((1/3) * acos(q_ / p^(1.5))) + b/3; % Eqn.10
  [~, n] = size(vb);

  if (n == 2)
    g3 = sqrt(2 * sqrt(d) - b); % Eqn.18
    g4 = sqrt(-2 * sqrt(d) - b); % Eqn.17
    lambda = (g3 + g4) / 2;
  else
    g1 = sqrt(u1 - b); % Eqn.15
    g2 = -2 * sqrt(u1^2 - 4*d); % check_me
    lambda = 0.5 * (g1 + sqrt(-u1 - b - g2)); % Eqn.14
  end

  H = K - lambda * eye(4); % Eqn.19
  qk = zeros(4, 1);
  q = zeros(4, 1);

  % k = i
  for i = 1 : 4
    H_ki = H;
    H_ki(i,:) = []; % Delete row
    H_ki(:,i) = []; % Delete column
    qk(i) = det(H_ki);
  end

  % Farthest-from-zero element
  [~,k] = max(abs(qk));

  % Optimal quaternion, Eqn.20
  for i = 1 : 4
    H_ki = H;
    H_ki(k,:) = [];
    H_ki(:,i) = [];
    q(i) = (-1)^(k+i) * det(H_ki);
  end

  q = q / norm(q);
  C = quaternion_to_dcm(q);
end
