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