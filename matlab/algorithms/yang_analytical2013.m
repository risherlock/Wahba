function [C_opt, q_opt] = yang_analytical2015(v_b, v_i, w)
%%% Yang's analytical solution to Wahba's problem
%
% Inputs:
%   v_b = Vectors in the spacecraft body frame i.e. measurements (3xn).
%   v_i = Corresponding vectors known in the inertial frame (3xn).
%
% Outputs:
%   q_opt = Optimal quaternion that transforms vectors in inertial frame to 
%           vectors in body frame ([q0,q1,q2,q3]').
%   C_opt = Optimal rotation matrix that transforms vectors in inertial 
%           frame to vectors in body frame.
%
% Reference:
%   [1] Yaguang Yang and Zhiqiang Zhou - An Analytic Solution to Wahba's 
%       Problem (2013)
%
% Rishav (2020-12-20)

% Eqn.3
B = (v_b.*repmat(w,[1,3])')*v_i';
S = B + B';
z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
K = [B + B'- eye(3)*trace(B), z; z', trace(B)]; 

adj_S = det(S)*inv(S);
adj_K = det(K)*inv(K);
adj_B = det(B)*inv(B);

% Characteristic quatric equation of K-matrix: (Eqn.6)
% p(x) = x^4 + a^3 + b*x^3 + c*x + d = 0
a = 0;
b = -2*(trace(B))^2 + trace(adj_S) - z'*z;
c = -trace(adj_K);
d = det(K);

% Analytical solution Eqn.6
[~,lambda] = quartic(a,b,c,d);

% Determine optimal DCM (Markley1993)
norm_B = norm(B, 'fro');
kappa = 0.5*(lambda^2 - norm_B^2);
zeta = kappa*lambda - det(B);

C_opt = ((kappa + norm_B^2)*B + lambda*adj_B' - B*(B')*B)/zeta;
q_opt = dcm_to_quaternion(C_opt);
end

function [x,x_max] = quartic(a,b,c,d)
%%% Computes roots of quatric equation analytically 
%
% Inputs:
%   a,b,c,d = Coefficient of the quartic equation: 
%               p(x) = x4 + ax3 + bx2 + cx + d = 0.
%
% Outputs:
%   x    = Roots of the quartic equation, [x1, x2, x3, x4]'.
%   xmax = Maximum real root of the quartic equation.
%
% References:
%   Yaguang Yang, Zhiqiang Zhou - An Analytic Solution to Wahba's Problem (2013)
%
% Rishav (2020/12/21)

% Tolerance for Eqn.13c
tolerance = 10e-3;

% Eqn.9
p = a*c - b^2/3 - 4*d;
q = a*b*c/3 - (a^2)*d - (2/27)*b^3 - c^2 + (8/3)*b*d;
w1 =  - 1 + 1i*sqrt(3)/2;
w2 =  - 1 - 1i*sqrt(3)/2;

% Cardano's formula (Eqn.10)
y1 = (- q/2 + sqrt((q/2)^2 + (p/3)^3))^(1/3) ...
   + (- q/2 - sqrt((q/2)^2 + (p/3)^3))^(1/3);
y2 = w1*(- q/2 + sqrt((q/2)^2 + (p/3)^3))^(1/3) ...
   + w2*(- q/2 - sqrt((q/2)^2 + (p/3)^3))^(1/3);
y3 = w2*(- q/2 + sqrt((q/2)^2 + (p/3)^3))^(1/3) ...
   + w1*(- q/2 - sqrt((q/2)^2 + (p/3)^3))^(1/3);

% Eqn.12 holds for real y
if imag(y1) == 0
    y = y1;
elseif imag(y2) == 0
    y = y2;
elseif imag(y3) == 0
    y = y3;
end

% Eqn.12
g1 = + sqrt(y - (2/3)*b);
g2 = - sqrt(y - (2/3)*b);
h1 = (y1 + b/3 + sqrt((y1 + b/3)^2 - 4*d))/2;
h2 = (y1 + b/3 - sqrt((y1 + b/3)^2 - 4*d))/2;

% Check if g1*h2 + g2*h1 = c holds.
% If not swap the variables
if abs(g1*h2+g2*h1-c) > tolerance
    h  = h1; 
    h1 = h2; 
    h2 = h;
end

% Roots of quartic (Eqn.14)
x1 = 0.5*(-g1 + sqrt(g1^2 - 4*h1));
x2 = 0.5*(-g1 - sqrt(g1^2 - 4*h1));
x3 = 0.5*(-g2 + sqrt(g2^2 - 4*h2));
x4 = 0.5*(-g2 - sqrt(g2^2 - 4*h2));

x  = [real(x1), real(x2), real(x3), real(x4)];
x_max = max(x);
end
