function [ RM ] = QtoRM( Q )
%QTORM Takes in quaternion versor Q and returns Rotation Matrix RM
%   Quaternion should be in the form [qr, qi, qj, qk] where qi, qj, and qk
%   represent the scalars dictating the vector around which we will rotate
%   the
%   NOTE - Limit use as much as possible. Degeneration is inherant in
%   conversion as the sin of the rotational angle approaches zero or as the
%   quaternion approaches the identity quanternion. Care has been taken in
%   the implementation of this program to limit this loss but loss will
%   inevitably occur.
%
% 2025-06-26 Ty Hedrick added check for identity quaternion

% check for identity quaternion
if Q == [1 0 0 0]
    RM = eye(3);
else

    a = 2*acos(Q(1));
    s = sin(a);
    c = cos(a);
    v = (1/(sin(a/2)))*Q(2:4);
    q = 1-c;

    RM = [c + (v(1)^2)*q, v(1)*v(2)*q-v(3)*s, v(1)*v(3)*q+v(2)*s;...
        v(2)*v(1)*q + v(3)*s, c + (v(2)^2)*q, v(2)*v(3)*q- v(1)*s;...
        v(3)*v(1)*q-v(2)*s, v(3)*v(2)*q+ v(1)*s, c+(v(3)^2)*q];
end

end