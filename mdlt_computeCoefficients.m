function [coefs, resid] = mdlt_computeCoefficients(xyz, imageuv)
%
%function [coefs,resid] = mdlt_computeCoefficients(xyz,imageuv)
%
% Inputs:
% xyz: [3,n] matrix of computed 3D positions
% imageuv: [2,n] matrix of corresponding image points
%
% REFERENCES: 1. The nonlinear constraint on the orthogonality that is used
% to solve for one of the coefficients in terms of the other 10 was taken
% from http://www.kwon3d.com/theory/dlt/mdlt.html
%
% 2. This paper shows the particular mathematical linearization technique
% for solving non-linear nature of equations due to adding non-linear
% constraints.
%
%	Miller N. R., Shapiro R., and McLaughlin T. M. A Technique for
%	Obtaining Spatial Kinematic Parameters of Segments of Biomechanical
%	Systems from Cinematographic Data. J. Biomech, 1980, v.13, pp535-547
%
% Ty Hedrick and Evan Bluhm

% check for any NaN rows (missing data) in the xyz or imageuv
ndx=find(sum(isnan([xyz,imageuv]),2)>0);

% remove any missing data rows
xyz(ndx,:)=[];
imageuv(ndx,:)=[];

% Compute initial linear least-squares estimate
M=zeros(size(xyz,1)*2,11);
for i=1:size(xyz,1)
  M(2*i-1,1:3)=xyz(i,1:3);
  M(2*i,5:7)=xyz(i,1:3);
  M(2*i-1,4)=1;
  M(2*i,8)=1;
  M(2*i-1,9:11)=xyz(i,1:3).*-imageuv(i,1);
  M(2*i,9:11)=xyz(i,1:3).*-imageuv(i,2);
end

% re-arrange the imageuv array for the linear solution
imageuvF=reshape(flipud(rot90(imageuv)),numel(imageuv),1);

% get the linear solution to the 11 parameters
coefs=M\imageuvF;

% do some iteration to ensure orthogonality: the hard part Using the
% orthogonality contstraint, estimate C1 The basic constraint:
% (C1*C5+C2*C6+C3*C7)*(C9^2+C10^2+C11^2)=(C1*C9+C2*C10+C3*C11) ...
%    *(C5*C9+C6*C10+C7*C11)
% Following (as did Pribanic) the linear least squares technique used in
% Miller (1980), solve for C2, ... C11 using this estimate.  Iterate.
maxIterations=30;
iter=1;
dC=999;
while iter<=maxIterations && dC>10^-7
  
  newCoefs = coefs;
  % Solve for C1 in terms of C2..C11
  newCoefs(1)= -((coefs(2)*coefs(10) + coefs(3)*coefs(11))*(coefs(5)*coefs(9) + ...
    coefs(6)*coefs(10) + coefs(7)*coefs(11)) - (coefs(2)*coefs(6) + ...
    coefs(3)*coefs(7))*(coefs(9)^2 + coefs(10)^2 + coefs(11)^2))/ ...
    (coefs(9)*(coefs(5)*coefs(9) + coefs(6)*coefs(10) + coefs(7)*coefs(11)) ...
    - coefs(5)*(coefs(9)^2 + coefs(10)^2 + coefs(11)^2));
  
  % Use the basic DLT equations
  % u = (C1*X + C2*Y + C3*Z) / (L9*X + L10*Y + L11*Z + 1)
  % v = (C5*X + C6*Y + C7*Z) / (L9*X + L10*Y + L11*Z + 1)
  %
  % and the new value for C1 to solve for C2, ..., C11
  
  % syms C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 X Y Z real
  % constraint = (C1*C5+C2*C6+C3*C7)*(C9^2+C10^2+C11^2)-(C1*C9+C2*C10+C3*C11) ...
  %    *(C5*C9+C6*C10+C7*C11);
  % C1=solve(constraint,C1);
  % u(X,Y,Z) = (C1*X + C2*Y + C3*Z +C4) / (C9*X + C10*Y + C11*Z + 1);
  % v(X,Y,Z) = (C5*X + C6*Y + C7*Z +C8) / (C9*X + C10*Y + C11*Z + 1);
  
  % Set up expression (9) in Miller with DLT formulation instead of explicit
  % camera extrinsics
  
  Zmat = zeros(2*size(imageuv,1),10);
  Fmat = zeros(2*size(imageuv,1),1);
  C2=newCoefs(2); C3=newCoefs(3); C4=newCoefs(4); C5=newCoefs(5);
  C6=newCoefs(6); C7=newCoefs(7); C8=newCoefs(8); C9=newCoefs(9);
  C10=newCoefs(10); C11=newCoefs(11);
  for i=1:size(imageuv,1)
    X = xyz(i,1); Y = xyz(i,2); Z = xyz(i,3);
    Zmat(i*2-1,1) = (Y - (X*(C10*(C5*C9 + C6*C10 + C7*C11) - C6*(C9^2 + ...
      C10^2 + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 ...
      + C11^2)))/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2-1,2) = (Z - (X*(C11*(C5*C9 + C6*C10 + C7*C11) - C7*(C9^2 + ...
      C10^2 + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 ...
      + C11^2)))/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2-1,3) = 1/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2-1,4) = -((C9*X*(C2*C10 + C3*C11))/(C9*(C5*C9 + C6*C10 + ...
      C7*C11) - C5*(C9^2 + C10^2 + C11^2)) + (X*(C10^2 + C11^2)*((C2*C10 ...
    + C3*C11)*(C5*C9 + C6*C10 + C7*C11) - (C2*C6 + C3*C7)*(C9^2 + C10^2 ...
    + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + ...
    C11^2))^2)/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2-1,5) = -((X*(C10*(C2*C10 + C3*C11) - C2*(C9^2 + C10^2 + ...
      C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + ...
      C11^2)) - (C9*C10*X*((C2*C10 + C3*C11)*(C5*C9 + C6*C10 + C7*C11) ...
      - (C2*C6 + C3*C7)*(C9^2 + C10^2 + C11^2)))/(C9*(C5*C9 + C6*C10 + ...
      C7*C11) - C5*(C9^2 + C10^2 + C11^2))^2)/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2-1,6) = -((X*(C11*(C2*C10 + C3*C11) - C3*(C9^2 + C10^2 + ...
      C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + ...
      C11^2)) - (C9*C11*X*((C2*C10 + C3*C11)*(C5*C9 + C6*C10 + C7*C11) ...
      - (C2*C6 + C3*C7)*(C9^2 + C10^2 + C11^2)))/(C9*(C5*C9 + C6*C10 ...
      + C7*C11) - C5*(C9^2 + C10^2 + C11^2))^2)/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2-1,7) = 0;
    Zmat(i*2-1,8) = ((X*(2*C9*(C2*C6 + C3*C7) - C5*(C2*C10 + C3*C11)))...
      /(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + C11^2)) + ...
      (X*(C6*C10 + C7*C11)*((C2*C10 + C3*C11)*(C5*C9 + C6*C10 + C7*C11) ...
      - (C2*C6 + C3*C7)*(C9^2 + C10^2 + C11^2)))/(C9*(C5*C9 + C6*C10 + ...
      C7*C11) - C5*(C9^2 + C10^2 + C11^2))^2)/(C9*X + C10*Y + C11*Z + 1)...
      - (X*(C4 + C2*Y + C3*Z - (X*((C2*C10 + C3*C11)*(C5*C9 + C6*C10 + ...
      C7*C11) - (C2*C6 + C3*C7)*(C9^2 + C10^2 + C11^2)))/(C9*(C5*C9 + ...
      C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + C11^2))))/(C9*X + C10*Y + ...
      C11*Z + 1)^2;
    Zmat(i*2-1,9) = - ((X*(C2*(C5*C9 + C6*C10 + C7*C11) - 2*C10*(C2*C6 ...
      + C3*C7) + C6*(C2*C10 + C3*C11)))/(C9*(C5*C9 + C6*C10 + C7*C11) ...
      - C5*(C9^2 + C10^2 + C11^2)) + (X*(2*C5*C10 - C6*C9)*((C2*C10 + ...
      C3*C11)*(C5*C9 + C6*C10 + C7*C11) - (C2*C6 + C3*C7)*(C9^2 + C10^2 ...
      + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + ...
      C11^2))^2)/(C9*X + C10*Y + C11*Z + 1) - (Y*(C4 + C2*Y + C3*Z -...
      (X*((C2*C10 + C3*C11)*(C5*C9 + C6*C10 + C7*C11) - (C2*C6 + ...
      C3*C7)*(C9^2 + C10^2 + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) ...
      - C5*(C9^2 + C10^2 + C11^2))))/(C9*X + C10*Y + C11*Z + 1)^2;
    Zmat(i*2-1,10) = - ((X*(C3*(C5*C9 + C6*C10 + C7*C11) - 2*C11*(C2*C6 ...
      + C3*C7) + C7*(C2*C10 + C3*C11)))/(C9*(C5*C9 + C6*C10 + C7*C11) ...
      - C5*(C9^2 + C10^2 + C11^2)) + (X*(2*C5*C11 - C7*C9)*((C2*C10 + ...
      C3*C11)*(C5*C9 + C6*C10 + C7*C11) - (C2*C6 + C3*C7)*(C9^2 + C10^2 ...
      + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + ...
      C11^2))^2)/(C9*X + C10*Y + C11*Z + 1) - (Z*(C4 + C2*Y + C3*Z - ...
      (X*((C2*C10 + C3*C11)*(C5*C9 + C6*C10 + C7*C11) - (C2*C6 + ...
      C3*C7)*(C9^2 + C10^2 + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) ...
      - C5*(C9^2 + C10^2 + C11^2))))/(C9*X + C10*Y + C11*Z + 1)^2;
    Zmat(i*2,1) = 0;
    Zmat(i*2,2) = 0;
    Zmat(i*2,3) = 0;
    Zmat(i*2,4) = X/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2,5) = Y/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2,6) = Z/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2,7) = 1/(C9*X + C10*Y + C11*Z + 1);
    Zmat(i*2,8) = -(X*(C8 + C5*X + C6*Y + C7*Z))/(C9*X + C10*Y + C11*Z + 1)^2;
    Zmat(i*2,9) = -(Y*(C8 + C5*X + C6*Y + C7*Z))/(C9*X + C10*Y + C11*Z + 1)^2;
    Zmat(i*2,10) = -(Z*(C8 + C5*X + C6*Y + C7*Z))/(C9*X + C10*Y + C11*Z + 1)^2;
    Fmat(i*2-1) = imageuv(i,1)-((C4 + C2*Y + C3*Z - (X*((C2*C10 + ...
      C3*C11)*(C5*C9 + C6*C10 + C7*C11) - (C2*C6 + C3*C7)*(C9^2 + C10^2 ...
      + C11^2)))/(C9*(C5*C9 + C6*C10 + C7*C11) - C5*(C9^2 + C10^2 + ...
      C11^2)))/(C9*X + C10*Y + C11*Z + 1));
    Fmat(i*2) = imageuv(i,2)-(C5*X + C6*Y + C7*Z +C8) / (C9*X + C10*Y + C11*Z + 1);
  end
  % Find least-squares solution
  delC=(Zmat'*Zmat)\(Zmat'*Fmat);
  coefs(2:11) = newCoefs(2:11)+delC;
  coefs(1)= newCoefs(1);
  % Check for convergence
  dC = norm(delC./newCoefs(2:11));
  iter=iter+1;
  
end

[uv]=dlt_inverse(coefs,xyz);
resid=(sum(sum((uv-imageuv).^2))./numel(imageuv))^0.5;
return
