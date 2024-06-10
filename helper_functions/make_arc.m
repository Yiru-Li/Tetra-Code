function p = make_arc(Tnodes, Tedges, A, B, C, extend)
arguments
    Tnodes (:, 3)
    Tedges (:, 2)
    A (1, 3)
    B (1, 3)
    C (1, 3)
    extend = false
end
arc_center = dot(B-C, A-C)/norm(A-C)*(A-C)/norm(A-C)+C;
xnew = A-arc_center;
xnew = xnew/norm(xnew);
ynew = B-arc_center;
ynew = ynew/norm(ynew);
znew = cross(xnew, ynew);
transform = [xnew; ynew; znew];
% search edges intersecting xy-plane
transformed_nodes = (Tnodes-arc_center)/transform;
intersecting_edges = Tedges(xor(transformed_nodes(Tedges(:, 1), 3)>0, transformed_nodes(Tedges(:, 2), 3)>0), :);
above = transformed_nodes(intersecting_edges(:, 1), :);
below = transformed_nodes(intersecting_edges(:, 2), :);
delta = above-below;
delta_z = above(:, 3)./delta(:, 3);
arc = above-(delta.*delta_z);
A_transform = (A-arc_center)/transform;
C_transform = (C-arc_center)/transform;
arc = [A_transform; arc; C_transform];
if extend
    arc = arc(:, 1:2);
    [~, p] = freeBoundary(delaunayTriangulation(arc(arc(:, 2)*3-arc(:, 1)>-100, :)));
    p = [p zeros(height(p), 1)];
    ID = dsearchn(p, A_transform);
    p = [A; p([ID+1:end 1:ID-1], :)*transform+arc_center];
else
    arc = arc(arc(:, 2)>0, :);
    [~, p] = freeBoundary(delaunayTriangulation(arc(:, 1:2)));
    p = [p zeros(height(p), 1)];
    IDs = dsearchn(p, [A_transform; C_transform]);
    if IDs(1)>IDs(2)
        p = p([IDs(1):end 1:IDs(2)], :);
    else
        p = p(IDs(1):IDs(2), :);
    end
    p = [A; p(2:end-1, :)*transform+arc_center; C];
end
p = unique(p,'rows', 'stable');
end