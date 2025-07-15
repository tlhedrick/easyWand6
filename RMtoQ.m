function [ Q ] = RMtoQ( RM )
%QTORM Takes in Rotation Matrix RM and returns quaternion versor Q.
%   RM should be in the form of a 3x3 rotational matrix of the form 
%   R =[ux, vx, wx]
%      [uy, vy, wy]
%      [uz, vz, wz]
%   NOTE - Limit use as much as possible. Degeneration is inherant in 
%   conversion as the sin of the rotational angle approaches zero or as the
%   quaternion approaches the identity quanternion. Care has been taken in
%   the implementation of this program to limit this loss but loss will
%   inevitably occur.

    a = (1/2)*sqrt(1 + RM(1,1) + RM(2,2) + RM(3,3));
    b = (1/2)*sqrt(1 + RM(1,1) - RM(2,2) - RM(3,3));
    
    if a<b
        a = (1/(4*b))*(RM(3,2) - RM(2,3));
        c = (1/(4*b))*(RM(1,2) + RM(2,1));
        d = (1/(4*b))*(RM(1,3) + RM(3,1));
        
    else
        b = (1/(4*a))*(RM(3,2) - RM(2,3));
        c = (1/(4*a))*(RM(1,3) - RM(3,1));
        d = (1/(4*a))*(RM(2,1) - RM(1,2));
    end
    
    Q = [a,b,c,d];
    

    
    
end