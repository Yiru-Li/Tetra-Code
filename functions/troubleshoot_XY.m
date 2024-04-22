for h = 95676%1:height(tmp)
    X = intersect(HC, [0 0; tmp(h, 1:2)*abs(max(200./tmp(h, 1:2)))]);
    X = [X(all(X~=0, 2), :) 0]; % intersection with head circumference
    seg = find(all(HC_vertices==X(1:2), 2)); % in case intersection lands on a vertex
    if isempty(seg)
        [~, seg] = max(sum(diff(sign(HC_vertices-X(1:2)))~=0, 2));
        X_percent(h) = sum(sqrt(sum(([HC_vertices(2:seg, :); X(1:2)]-HC_vertices(1:seg, :)).^2, 2)))/X_total;
    else
        X_percent(h) = sum(sqrt(sum(([HC_vertices(2:seg, :)]-HC_vertices(1:seg-1, :)).^2, 2)))/X_total;
    end
    p = make_arc(Tnodes_new, Tedges, Cz_new, [X(1:2) Cz_new(3)], X);
    Y_total = sum(sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2)));
    if sqrt(sum((tmp(h, :)-X).^2))<1 % point close to head circumference
        Y_percent(h) = sum(sqrt(sum((p(1:end-1, :)-[p(2:end-1, :); tmp(h, :)]).^2, 2)))/Y_total;
    elseif tmp(h, 3)<0
        p = make_arc(Tnodes_new, Tedges(all(ismember(Tedges, find(Tnodes_new(:, 3)>tmp(h, 3))), 2), :), Cz_new, [X(1:2)/2 Cz_new(3)], tmp(h, :));
        Y_percent(h) = sum(sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2)))/Y_total;
    else
        [~, seg] = max(sum(diff(sign(p-tmp(h, :)))~=0, 2));
        if isempty(seg)
            [dist, id] = min(sqrt(sum(([Cz_new; X]-tmp(h , :)).^2, 2)));
            if id==1
                Y_percent(h) = dist/Y_total;
            else
                Y_percent(h) = sum(sqrt(sum((p(1:end, :)-[p(2:end, :); tmp(h, :)]).^2, 2)))/Y_total;
            end
        else
            assert(length(seg)<=1, 'multiple matches');
            Y_percent(h) = sum(sqrt(sum(([p(2:seg, :); tmp(h, :)]-p(1:seg, :)).^2, 2)))/Y_total;
        end
    end
    patch('Faces',scalp.triangles,'Vertices',Tnodes_new);
    hold on
    plot3(p(:, 1), p(:, 2), p(:, 3))
    scatter3(tmp(h, 1), tmp(h, 2), tmp(h, 3))
    hold off
end
%% Function definitions
function [refs, arc_lengths] = extend_arc(p, ref_length, ref_points, total_points)
dists = sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2));
cum_dists = [0; cumsum(dists)];
arc_lengths = linspace(0, ref_length/(ref_points-1)*(total_points-1), total_points);
refs = zeros(total_points, 3);
refs(1, :) = p(1, :);
for n = 2:length(arc_lengths)
    startID = find(cum_dists<arc_lengths(n), 1, 'last');
    A = p(startID, :);
    B = p(startID+1, :);
    refs(n, :) = A+(B-A)*(arc_lengths(n)-cum_dists(startID))/dists(startID);
end
end

function [refs, arc_lengths] = divide_arc(p, num_arc_refs)
dists = sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2));
cum_dists = [0; cumsum(dists)];
arc_lengths = linspace(0, cum_dists(end), num_arc_refs);
refs = zeros(num_arc_refs, 3);
refs(1, :) = p(1, :);
refs(end, :) = p(end, :);
for n = 2:length(arc_lengths)-1
    startID = find(cum_dists<arc_lengths(n), 1, 'last');
    A = p(startID, :);
    B = p(startID+1, :);
    refs(n, :) = A+(B-A)*(arc_lengths(n)-cum_dists(startID))/dists(startID);
end
end

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
end