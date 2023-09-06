function q = dcm_to_quaternion(C)
%%% Quaternion representated by Direction Cosine Matrix.
%
% Reference:
%   Malcolm D. Shuster - A Survey of Attitude Representations (1993)
%
% Rishav (2020-12-21)

% Tolerance of q0
tolerance = 1e-5;

% q0 ~ 0 is numerically insignificant in square root (Shuster1993)
q0 = sqrt(0.25*(trace(C)+1));

if (q0 > tolerance)
    q0x4_inv = 0.25/q0; 
    q1 = q0x4_inv*(C(2,3) - C(3,2));
    q2 = q0x4_inv*(C(3,1) - C(1,3));
    q3 = q0x4_inv*(C(1,2) - C(2,1));
else
    q1 = sqrt(0.25*(1 + C(1,1) - C(2,2) - C(3,3)));
    q1x4_inv = 0.25/q1;
    q2 = q1x4_inv * (C(1,2) + C(2,1));
    q3 = q1x4_inv * (C(1,3) + C(3,1));
    q0 = q1x4_inv * (C(2,3) - C(3,2));
end

q = [q1, q2, q3, q0]';
end
