function [y]=binarycombinations(n)

% function [y]=binarycombinations(n)
%
% Create a matrix of 2^n rows and n columns with a unique binary state in
% each row.  The first row is always all zeros and the final row all ones.
% Works through the possibilities in order, with all possibilities
% including only one "1" occuring before any possibilities with two "1"s
% and so on.
%
% Ty Hedrick

num_hyp=2^n;			% generate the set of hypotheses
y=zeros(n,num_hyp);		% all possible bit combos
for index=1:n,
  y(index,:)=(-1).^ceil([1:num_hyp]/(2^(index-1)));
end

% change -1s to 0s and rotate
idx=find(y==-1);
y(idx)=0;
y=y';

% sort
y(:,end+1)=y*ones(n,1);
y=sortrows(y,n+1);
y(:,end)=[];