function [C, q] = triad1964(b1,b2,r1,r2)
%%% Tri-axial Attitude Determination (TRIAD)
%
% Inputs:
%   b1 = First unit vector in spacecraft body frame (1x3)
%   b2 = Second unit vector in spacecraft body frame (1x3)
%   r1 & r2 = Vectors corresponding to b1 & b2 in the reference frame
% 
% Output:
%   C = Rotation matrix that transforms vectors in inertial frame to the 
%       vectors in body frame.
%   q = Quaternion that transforms vectors in inertial frame to the vectors 
%        in body frame.
%
% Note:
%   b1 is the one that is much more accurately determined than b2, so the 
%   estimate satisfies b1 = Q*r1 exactly, but b2 = Q*r2 only approximately.
%
% Reference: 
%   [1] Harold D. Black - A Passive System for Determining the Attitude of 
%       a Satellite (1964)
%  
% Rishav (2020-11-4)

% Triads vectors {v1,v2,v3} using r1 & r2
% v1 = r1; 
% v3 = cross(r1,v2);
v2 = cross(r1,r2)/norm(cross(r1,r2));

% Triads vectors {w1,w2,w2}' using r1&r2
% w1 = b1; 
% w3 = cross(b1,w2);
w2 = cross(b1,b2)/norm(cross(b1,b2));

C = b1*r1' + cross(b1,w2)*cross(r1,v2)' + w2*v2';
q = dcm_to_quaternion(C);
end
