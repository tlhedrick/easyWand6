function  [xyz, reprojErrors] = triangulate_v4(Ks, Rs, ts, uv)

% function [xyz, reprojErrors] = triangulate_v4(Ks, Rs, ts, uv)
%
% Inputs:
%   Ks          - 1xN cell array of 3x3 intrinsics matrices
%   Rs          - 1xN cell array of 3x3 rotation matrices
%   ts          - 1xN cell array of 3x1 translation vectors
%   uv      - MxNx2 array of image points: M points, N cameras, (u, v)
%   valid_mask  - MxN logical array indicating valid observations (1: valid, 0: occluded)
%
% Outputs:
%   xyz          - Mx3 matrix of estimated 3D points
%   reprojErrors - Mx1 vector of average reprojection error per point
%
% Linear multi-camera reprojection code developed from a ChatGPT baseline.
% 
% Example K for formatting info:
%  K = [500 0 320; 0 500 240; 0 0 1]; --> Not Mathworks style
%
% Ty Hedrick, 2025-06-24

% setup outputs
[M, N, ~] = size(uv); % M points, N/2 cameras
N=N/2; % N cameras
xyz = zeros(M, 3);
reprojErrors = zeros(M, 1);

% Precompute projection matrices
Ps = cell(1, N);
for i = 1:N
    Ps{i} = Ks{i} * [Rs{i}, ts{i}]; % testing -ts{i} to match Mathworks --> totally blows things up
    % cameraMatrices{id} = (intrinsics(i).K * [R, -R*t'])';
end

for ptIdx = 1:M % for each point
    % Collect valid observations for this point
    A = [];
    valid_cam_count = 0;  % Track how many valid cameras there are for this point
    for camIdx = 1:N
        if isfinite(uv(ptIdx,camIdx*2)) % check for good data from this camera
            u = uv(ptIdx, camIdx*2-1);
            v = uv(ptIdx, camIdx*2);
            P = Ps{camIdx};

            A = [A;
                u * P(3,:) - P(1,:);
                v * P(3,:) - P(2,:)];

            valid_cam_count = valid_cam_count + 1;
        end
    end

    % Proceed only if there are valid cameras
    if valid_cam_count >= 2
        % Solve using SVD
        [~, ~, V] = svd(A);
        xyz_hom = V(:, end);
        xyz_hom = xyz_hom / xyz_hom(4); % change to -xyz_hom(4) on 2025-06-30 to get a right-hand result??
        xyz(ptIdx, 1:3) = xyz_hom(1:3)';

        % Reproject and compute error
        totalError = 0;
        for camIdx = 1:N
            if isfinite(uv(ptIdx,camIdx*2))  % Only reproject for valid cameras
                P = Ps{camIdx};
                xyz_proj = P * [xyz_hom(1:3); 1]; % re-homogenize
                xyz_proj = xyz_proj ./ xyz_proj(3);
                observed = squeeze(uv(ptIdx, camIdx, :));
                observed = squeeze(uv(ptIdx,camIdx*2-1:camIdx*2,:));
                error = norm(observed - xyz_proj(1:2)');
                totalError = totalError + error;
            end
        end

        reprojErrors(ptIdx) = totalError / valid_cam_count;
    else
        % If less than 2 valid observations, set the point to NaN
        xyz(ptIdx, :) = NaN;
        reprojErrors(ptIdx) = NaN;
    end
end
end

