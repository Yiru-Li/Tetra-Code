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