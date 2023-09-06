function [C] = quaternion_to_dcm(q)
%%% DCM in terms of quaternions
% Maps a vector from inertial to body frame: V_b = C * V_i
%   
% Inputs:
%   q = [q1, q2, q3,  q0]', a quaternion vector   
%
% Output:
%   C = Direction cosine matrix
%
% Reference:
%   Schaub - Analytical Mechanics of Space Systems (2nd ed.)(Pg.89)
%
% Rishav (2020-6-23)

q0 = q(4);
q1 = q(1);
q2 = q(2);
q3 = q(3);

C(1,1) = q0^2 + q1^2 - q2^2 - q3^2;
C(1,2) = 2*(q1*q2 + q0*q3);
C(1,3) = 2*(q1*q3 - q0*q2);
C(2,1) = 2*(q1*q2 - q0*q3);
C(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3;
C(2,3) = 2*(q2*q3 + q0*q1);
C(3,1) = 2*(q1*q3 + q0*q2);
C(3,2) = 2*(q2*q3 - q0*q1);
C(3,3) = q0^2 - q1^2 - q2^2 + q3^2;
end
