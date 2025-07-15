function [R, t] = twoCamCal_v4(E, K1, K2, x1, x2)
%
% function [R, t] = twoCamCal_v4(E, K1, K2, x1, x2)
%
% Compute R and t from Essential matrix, intrinsics, and a set of shared
% observations
%
% Inputs:
%   E  - Essential matrix (3x3)
%   K1, K2 - Camera intrinsics (3x3 each)
%   x1, x2 - Nx2 corresponding image points in image 1 and 2
% Outputs:
%   R - Relative rotation from camera 1 to 2 (3x3)
%   t - Relative translation from camera 1 to 2 (3x1, unit vector)

    % Normalize image points (convert to camera coordinates)
    x1_norm = (K1 \ [x1, ones(size(x1,1),1)]')';
    x2_norm = (K2 \ [x2, ones(size(x2,1),1)]')';

    % SVD decomposition of E
    [U, ~, V] = svd(E);

    % Ensure proper rotation matrices
    if det(U) < 0, U = -U; end
    if det(V) < 0, V = -V; end

    % W matrix as per Hartley & Zisserman
    W = [0 -1 0; 1 0 0; 0 0 1];

    % Four possible decompositions
    R1 = U * W  * V';
    R2 = U * W' * V';
    t1 = U(:,3);
    t2 = -U(:,3);

    % Ensure rotations are proper
    if det(R1) < 0, R1 = -R1; end
    if det(R2) < 0, R2 = -R2; end

    % Test all 4 combinations
    candidates = {
        R1, t1;
        R1, t2;
        R2, t1;
        R2, t2
    };

    % P1 is canonical: [I | 0]
    P1 = [eye(3), zeros(3,1)];

    max_positive_depth = 0;
    for i = 1:4
        R_test = candidates{i, 1};
        t_test = candidates{i, 2};
        P2 = [R_test, t_test];

        % Triangulate points
        pts3D = linearTriangulation(P1, P2, x1_norm, x2_norm);

        % Cheirality check: points must be in front of both cameras
        num_positive = sum(pts3D(:,3) > 0 & (pts3D * R_test(3,:)' + t_test(3)) > 0);

        if num_positive > max_positive_depth
            max_positive_depth = num_positive;
            R = R_test;
            t = t_test;
        end
    end

    % Normalize t (optional: direction only)
    t = t / norm(t);
    % t = -t / norm(t); was this on 2025-06-30
end

function X = linearTriangulation(P1, P2, x1, x2)
% LINEARTRIANGULATION Triangulate points from two views
    N = size(x1,1);
    X = zeros(N,3);

    for i = 1:N
        A = [
            x1(i,1)*P1(3,:) - P1(1,:);
            x1(i,2)*P1(3,:) - P1(2,:);
            x2(i,1)*P2(3,:) - P2(1,:);
            x2(i,2)*P2(3,:) - P2(2,:);
        ];
        [~, ~, V] = svd(A);
        X_h = V(:,end);
        X(i,:) = (X_h(1:3) / X_h(4))';
    end
end
