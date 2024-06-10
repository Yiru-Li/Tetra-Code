clear all
addpath('./helper_functions/')
addpath('./functions_from_SimNIBS/')
lookup_initial_fit = readtable('Fiducials.csv');
m = mesh_load_gmsh4('MNI152.msh');
scalp = mesh_extract_regions(m, 'region_idx', 1005);
T = triangulation(double(scalp.triangles), scalp.nodes);
Tedges = edges(T);
az_steps = 400;
az_steps_NzIz = 361;
elev_steps = 399;
[tet_az,tet_elev] = meshgrid(1:az_steps, 1:elev_steps);
tet_az = reshape(tet_az, [], 1);
tet_elev = reshape(tet_elev, [], 1);
tet = array2table([tet_az tet_elev zeros(az_steps*elev_steps, 3)]);
tet.Properties.VariableNames = {'tet_az', 'tet_elev', 'MNI_x', 'MNI_y', 'MNI_z'};
lookup_initial_fit.Properties.VariableNames = {'Label', 'x', 'y', 'z', 'Name'};
%% finding the true location of Cz
Cz = scalp.nodes(scalp.nodes(:, 3)==max(scalp.nodes(:, 3)), :); % Cz estimate
LPA = lookup_initial_fit{strcmp(lookup_initial_fit.Name, 'LPA'), ["x", "y", "z"]}; % start point
RPA = lookup_initial_fit{strcmp(lookup_initial_fit.Name, 'RPA'), ["x", "y", "z"]}; % end point
Nz = lookup_initial_fit{strcmp(lookup_initial_fit.Name, 'Nz'), ["x", "y", "z"]}; % start point
Iz = lookup_initial_fit{strcmp(lookup_initial_fit.Name, 'Iz'), ["x", "y", "z"]}; % end point
Cz_eps = 1;
while Cz_eps>1e-4
    % arc fitting along Nz-Iz line
    NzIz = make_arc(scalp.nodes, Tedges, Nz, Cz, Iz);
    [NzIz, arc_lengths] = divide_arc(NzIz, 11);
    Cz = NzIz(6, :);
    % arc fitting along LPA-RPA line
    LPARPA = make_arc(scalp.nodes, Tedges, LPA, Cz, RPA);
    LPARPA = divide_arc(LPARPA, 11);
    Cz = LPARPA(6, :);
    Cz_eps = sum((Cz-NzIz(6, :)).^2);
end
NzIz_length = arc_lengths(end);
%% extension
NzIz_full = make_arc(scalp.nodes, Tedges, Nz, Cz, Iz, true);
NzIz_full = extend_arc(NzIz_full, NzIz_length, az_steps_NzIz, az_steps);
% mesh_show_surface(m, 'region_idx', [1002 1005], 'facealpha', 0.5)
% hold on
% plot3(NzIz_full(:, 1), NzIz_full(:, 2), NzIz_full(:, 3), 'LineWidth', 5)
% plot3(NzIz_full(:, 1), NzIz_full(:, 2), NzIz_full(:, 3), 'LineWidth', 5)
% view(90, 0)
%% initial fit
for az = 1:az_steps
    p = make_arc(scalp.nodes, Tedges, LPA, NzIz_full(az, :), RPA);
    refs = divide_arc(p, elev_steps);
    tet{tet.tet_az==az, {'MNI_x', 'MNI_y', 'MNI_z'}} = refs;
end
%% quantify scalp resolution
x = reshape(tet.MNI_x, elev_steps, az_steps);
y = reshape(tet.MNI_y, elev_steps, az_steps);
z = reshape(tet.MNI_z, elev_steps, az_steps);
figure
t = tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile
az_res = sqrt(diff(x, 1, 2).^2+diff(y, 1, 2).^2+diff(z, 1, 2).^2);
histogram(az_res, 0:0.05:1.5, 'Normalization', 'probability')
xline(median(az_res,"all"))
text(median(az_res,"all"), 0.1, 'Median ', 'HorizontalAlignment', 'right')
title('Nz-Iz step size')
disp(['Maximum distance between longitudinal nodes = ' num2str(max(az_res, [], 'all')) ' mm'])
nexttile
elev_res = sqrt(diff(x, 1, 1).^2+diff(y, 1, 1).^2+diff(z, 1, 1).^2);
histogram(elev_res, 0:0.05:1.5, 'Normalization', 'probability')
xline(median(elev_res,"all"))
text(median(elev_res,"all"), 0.2, 'Median ', 'HorizontalAlignment', 'right')
title('LPA-RPA step size')
disp(['Maximum distance between lateral nodes = ' num2str(max(elev_res, [], 'all')) ' mm'])
xlabel(t, 'Distance (mm)')
ylabel(t, 'Proportion')
nexttile
coverage = max(max(sqrt(az_res(1:end-1, :).^2+elev_res(:, 1:end-1).^2), ...
    sqrt(az_res(1:end-1, :).^2+elev_res(:, 2:end).^2)), ...
    max(sqrt(az_res(2:end, :).^2+elev_res(:, 1:end-1).^2), ...
    sqrt(az_res(2:end, :).^2+elev_res(:, 2:end).^2)))/2; % half the diagonal distance
histogram(coverage, 0:0.05:1.5, 'Normalization', 'probability')
xline(median(coverage,"all"))
text(median(coverage,"all"), 0.1, ' Median')
disp(['Maximum distance to nearest Tetra code coordinate = ' num2str(max(coverage, [], 'all')) ' mm'])
title('Maximum distance to nearest Tetra code coordinate')
set(gcf, 'Units', 'inches', 'Position', [0 0 6.5 4])
%% Direct projection to brain
% if isfile('extended_brain_projection.mat')
%     load('extended_brain_projection.mat');
% else
%     tmp = tet{:, {'MNI_x', 'MNI_y', 'MNI_z'}};
%     nhat = cross(brain.nodes(brain.triangles(:, 2), :)-brain.nodes(brain.triangles(:, 1), :), brain.nodes(brain.triangles(:, 3), :)-brain.nodes(brain.triangles(:, 1), :), 2);
%     nhat = nhat./sqrt(dot(nhat, nhat, 2));
%     e1 = cross(nhat, brain.nodes(brain.triangles(:, 2), :)-brain.nodes(brain.triangles(:, 1), :), 2);
%     e2 = cross(nhat, brain.nodes(brain.triangles(:, 3), :)-brain.nodes(brain.triangles(:, 2), :), 2);
%     e3 = cross(nhat, brain.nodes(brain.triangles(:, 1), :)-brain.nodes(brain.triangles(:, 3), :), 2);
%     projection = zeros(length(tmp), 3);
%     [N1,N2,N3] = surfnorm(x, y, z);
%     N = [reshape(-N1, [], 1) reshape(-N2, [], 1) reshape(-N3, [], 1)];
%     for n = 1:length(tmp)
%         ahat = N(n, :);
%         pcloud = brain.nodes-tmp(n, :);
%         p1 = pcloud(brain.triangles(:, 1), :);
%         p2 = pcloud(brain.triangles(:, 2), :);
%         p3 = pcloud(brain.triangles(:, 3), :);
%         d = mean([dot(nhat, p1, 2) dot(nhat, p2, 2) dot(nhat, p3, 2)], 2)./dot(repmat(ahat, length(nhat), 1), nhat, 2); % rounding error
%         p = d.*ahat;
%         intersections = (dot(e1, p-p1, 2)>=0) & (dot(e2, p-p2, 2)>=0) & (dot(e3, p-p3, 2)>=0);
%         mindis = min(abs(d(intersections)));
%         if ~isempty(mindis)
%             projection(n, :) = p((abs(d)==mindis) & intersections, :)+tmp(n, :);
%         else
%             projection(n, :) = nan(1, 3);
%         end
%     end
%     mesh_show_surface(brain)
%     hold on
%     scatter3(projection(:, 1), projection(:, 2), projection(:, 3), '.')
%     hold off
% end
% tet.Brain_x = projection(:, 1);
% tet.Brain_y = projection(:, 2);
% tet.Brain_z = projection(:, 3);
%% Balloon projection to brain
% median(sqrt(sum(((scalp.nodes(Tedges(:, 1), :)-scalp.nodes(Tedges(:, 2), :)).^2), 2))); % distance between nodes is about twice the Okamoto & Dan paper's
if isfile('extended_brain_projection.mat')
    load('extended_brain_projection.mat');
else
    brain = mesh_extract_regions(m, 'region_idx', 1002);
    tmp = tet{:, {'MNI_x', 'MNI_y', 'MNI_z'}};
    nhat = cross(brain.nodes(brain.triangles(:, 2), :)-brain.nodes(brain.triangles(:, 1), :), brain.nodes(brain.triangles(:, 3), :)-brain.nodes(brain.triangles(:, 1), :), 2);
    nhat = nhat./sqrt(dot(nhat, nhat, 2));
    e1 = cross(nhat, brain.nodes(brain.triangles(:, 2), :)-brain.nodes(brain.triangles(:, 1), :), 2);
    e2 = cross(nhat, brain.nodes(brain.triangles(:, 3), :)-brain.nodes(brain.triangles(:, 2), :), 2);
    e3 = cross(nhat, brain.nodes(brain.triangles(:, 1), :)-brain.nodes(brain.triangles(:, 3), :), 2);
    projection = zeros(length(tmp), 3);
    for n = 1:length(tmp)
        Idx = knnsearch(brain.nodes, tmp(n, :), 'K', 50);
        ahat = mean(brain.nodes(Idx, :))-tmp(n, :);
        ahat = ahat/norm(ahat);
        pcloud = brain.nodes-tmp(n, :);
        p1 = pcloud(brain.triangles(:, 1), :);
        p2 = pcloud(brain.triangles(:, 2), :);
        p3 = pcloud(brain.triangles(:, 3), :);
        d = mean([dot(nhat, p1, 2) dot(nhat, p2, 2) dot(nhat, p3, 2)], 2)./dot(repmat(ahat, length(nhat), 1), nhat, 2); % rounding error
        p = d.*ahat;
        intersections = (dot(e1, p-p1, 2)>=0) & (dot(e2, p-p2, 2)>=0) & (dot(e3, p-p3, 2)>=0);
        mindis = min(abs(d(intersections)));
        if ~isempty(mindis)
            projection(n, :) = p((abs(d)==mindis) & intersections, :)+tmp(n, :);
        else
            projection(n, :) = nan(1, 3);
        end
    end
    mesh_show_surface(brain)
    hold on
    scatter3(projection(:, 1), projection(:, 2), projection(:, 3), '.')
    hold off
end
tet.Brain_x = projection(:, 1);
tet.Brain_y = projection(:, 2);
tet.Brain_z = projection(:, 3);
%% brain coverage
mesh_show_surface(m, 'region_idx', 1002)
hold on
% mesh(x, y, z, 'FaceAlpha', 0)
scatter3(projection(:, 1), projection(:, 2), projection(:, 3), '.')
hold off
view(90, 0)
%% Create Tetra Codes
OLC_lookup = {'2' '3' '4' '5' '6' '7' '8' '9' 'C' 'F' 'G' 'H' 'J' 'M' 'P' 'Q' 'R' 'V' 'W' 'X'};
tet.tet_az = tet.tet_az-1;
tet.tet_elev = tet.tet_elev-1;
tet.OLC = [char(OLC_lookup(floor(tet.tet_az/20)+1).')...
    char(OLC_lookup(mod(tet.tet_az, 20)+1).')...
    char(OLC_lookup(floor(tet.tet_elev/20)+1).')...
    char(OLC_lookup(mod(tet.tet_elev, 20)+1).')];
%% Calculate SGP coordinates
tet.Pnz = tet.tet_az./(az_steps_NzIz-1);
tet.Pal = tet.tet_elev./(elev_steps-1);
%% calculate 10-10 electrodes and match to nearest Tetra Code
tet.EEG = repmat({'N/A'}, height(tet), 1);
tet.EEG(tet.Pal==0) = {'LPA'};
tet.EEG(tet.Pal==1) = {'RPA'};
EEG_coords = NzIz;
EEG_names = {'Nz' 'Fpz' 'AFz' 'Fz' 'FCz' 'Cz' 'CPz' 'Pz' 'POz' 'Oz' 'Iz'}.';
EEG_coords = [EEG_coords; LPARPA];
EEG_names = [EEG_names; {'LPA' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'RPA'}.'];

est = divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Fpz'), :), EEG_coords(strcmp(EEG_names, 'T7'), :), EEG_coords(strcmp(EEG_names, 'Oz'), :)), 5);
[EEG_names, ia] = unique([EEG_names; {'Fpz' 'Fp1' 'AF7' 'F7' 'FT7', 'T7'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Fpz'), :), est(2, :), EEG_coords(strcmp(EEG_names, 'T7'), :)), 6)];
EEG_coords = EEG_coords(ia, :);
[EEG_names, ia] = unique([EEG_names; {'T7' 'TP7' 'P7' 'PO7' 'O1' 'Oz'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'T7'), :), est(4, :), EEG_coords(strcmp(EEG_names, 'Oz'), :)), 6)];
EEG_coords = EEG_coords(ia, :);

est = divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Fpz'), :), EEG_coords(strcmp(EEG_names, 'T8'), :), EEG_coords(strcmp(EEG_names, 'Oz'), :)), 5);
[EEG_names, ia] = unique([EEG_names; {'Fpz' 'Fp2' 'AF8' 'F8' 'FT8', 'T8'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Fpz'), :), est(2, :), EEG_coords(strcmp(EEG_names, 'T8'), :)), 6)];
EEG_coords = EEG_coords(ia, :);
[EEG_names, ia] = unique([EEG_names; {'T8' 'TP8' 'P8' 'PO8' 'O2' 'Oz'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'T8'), :), est(4, :), EEG_coords(strcmp(EEG_names, 'Oz'), :)), 6)];
EEG_coords = EEG_coords(ia, :);

lobes_med = {'F' 'FC' 'CP' 'P'};
lobes_lat = {'F' 'FT' 'TP' 'P'};
for lobe = 1:4
    est = divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, [lobes_lat{lobe} '8']), :), EEG_coords(strcmp(EEG_names, [lobes_med{lobe} 'z']), :), EEG_coords(strcmp(EEG_names, [lobes_lat{lobe} '7']), :)), 5);
    [EEG_names, ia] = unique([EEG_names; cellstr([[lobes_lat{lobe}; repmat(lobes_med{lobe}, 4, 1)] ['8' '6' '4' '2' 'z'].'])], 'stable');
    EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, [lobes_lat{lobe} '8']), :), est(2, :), EEG_coords(strcmp(EEG_names, [lobes_med{lobe} 'z']), :)), 5)];
    EEG_coords = EEG_coords(ia, :);
    [EEG_names, ia] = unique([EEG_names; cellstr([[repmat(lobes_med{lobe}, 4, 1); lobes_lat{lobe}] ['z' '1' '3' '5' '7'].'])], 'stable');
    EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, [lobes_med{lobe} 'z']), :), est(4, :), EEG_coords(strcmp(EEG_names, [lobes_lat{lobe} '7']), :)), 5)];
    EEG_coords = EEG_coords(ia, :);
end

for lobe = {'AF' 'PO'}
    est = divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, [lobe{:} '8']), :), EEG_coords(strcmp(EEG_names, [lobe{:} 'z']), :), EEG_coords(strcmp(EEG_names, [lobe{:} '7']), :)), 5);
    [EEG_names, ia] = unique([EEG_names; cellstr([repmat(lobe{:}, 3, 1) ['8' '4' 'z'].'])], 'stable');
    EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, [lobe{:} '8']), :), est(2, :), EEG_coords(strcmp(EEG_names, [lobe{:} 'z']), :)), 3)];
    EEG_coords = EEG_coords(ia, :);
    [EEG_names, ia] = unique([EEG_names; cellstr([repmat(lobe{:}, 3, 1) ['z' '3' '7'].'])], 'stable');
    EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, [lobe{:} 'z']), :), est(4, :), EEG_coords(strcmp(EEG_names, [lobe{:} '7']), :)), 3)];
    EEG_coords = EEG_coords(ia, :);
end

est = divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Nz'), :), EEG_coords(strcmp(EEG_names, 'LPA'), :), EEG_coords(strcmp(EEG_names, 'Iz'), :)), 5);
[EEG_names, ia] = unique([EEG_names; {'Nz' 'N1' 'AF9' 'F9' 'FT9', 'LPA'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Nz'), :), est(2, :), EEG_coords(strcmp(EEG_names, 'LPA'), :)), 6)];
EEG_coords = EEG_coords(ia, :);
[EEG_names, ia] = unique([EEG_names; {'LPA' 'TP9' 'P9' 'PO9' 'I1' 'Iz'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'LPA'), :), est(4, :), EEG_coords(strcmp(EEG_names, 'Iz'), :)), 6)];
EEG_coords = EEG_coords(ia, :);

est = divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Nz'), :), EEG_coords(strcmp(EEG_names, 'RPA'), :), EEG_coords(strcmp(EEG_names, 'Iz'), :)), 5);
[EEG_names, ia] = unique([EEG_names; {'Nz' 'N2' 'AF10' 'F10' 'FT10', 'RPA'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'Nz'), :), est(2, :), EEG_coords(strcmp(EEG_names, 'RPA'), :)), 6)];
EEG_coords = EEG_coords(ia, :);
[EEG_names, ia] = unique([EEG_names; {'RPA' 'TP10' 'P10' 'PO10' 'I2' 'Iz'}.'], 'stable');
EEG_coords = [EEG_coords; divide_arc(make_arc(scalp.nodes, Tedges, EEG_coords(strcmp(EEG_names, 'RPA'), :), est(4, :), EEG_coords(strcmp(EEG_names, 'Iz'), :)), 6)];
EEG_coords = EEG_coords(ia, :);
[k, d] = dsearchn(tet{:, {'MNI_x', 'MNI_y', 'MNI_z'}}, EEG_coords);
tet.EEG(k(d<1)) = EEG_names(d<1);
mesh_show_surface(m, 'region_idx', [1002 1005], 'facealpha', 0.5)
hold on
mesh(x, y, z, 'FaceAlpha', 0, 'EdgeAlpha', 0.5)
plot3(EEG_coords(d<1&EEG_coords(:, 1)>0, 1), EEG_coords(d<1&EEG_coords(:, 1)>0, 2), EEG_coords(d<1&EEG_coords(:, 1)>0, 3), 'k.', 'MarkerSize', 18)
% text(EEG_coords(d<1&EEG_coords(:, 1)>0, 1)*1.2, EEG_coords(d<1&EEG_coords(:, 1)>0, 2)*1.1, EEG_coords(d<1&EEG_coords(:, 1)>0, 3)*1.1, EEG_names(d<1&EEG_coords(:, 1)>0), 'FontSize', 24)
hold off
view(90, 0)
%% Calculate X and Y
if ~exist('X_percent', 'var')
    Fpz = EEG_coords(strcmp(EEG_names, 'Fpz'), :);
    Oz = EEG_coords(strcmp(EEG_names, 'Oz'), :);
    sphere_center = dot(Cz-Fpz, Oz-Fpz)/norm(Oz-Fpz)*(Oz-Fpz)/norm(Oz-Fpz)+Fpz;
    xnew = (Fpz-sphere_center)/norm(Fpz-sphere_center);
    znew = (Cz-sphere_center)/norm(Cz-sphere_center);
    ynew = -cross(xnew, znew);
    transform = eye(3)/[xnew.' ynew.' znew.'];
    tmp = (transform*(tet{:, {'MNI_x', 'MNI_y', 'MNI_z'}}-sphere_center).').';
    Cz_new = (transform*(Cz-sphere_center).').';
    Fpz_new = (transform*(Fpz-sphere_center).').';
    Oz_new = (transform*(Oz-sphere_center).').';
    Tnodes_new = (transform*(scalp.nodes-sphere_center).').';
    X_total = 0;
    for B = [max(Tnodes_new(:, 2)) min(Tnodes_new(:, 2))]
        p = make_arc(Tnodes_new, Tedges, Fpz_new, [0 B 0], Oz_new);
        X_total = X_total + sum(sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2)));
    end
    X_percent = zeros(height(tmp), 1);
    Y_percent = zeros(height(tmp), 1);
    empty_list = [];
    for h = 1:height(tmp)
        if tmp(h, 3)<0
            p = make_arc(Tnodes_new, Tedges, Cz_new, [tmp(h, 1:2) 45], tmp(h, :));
        else
            l = sqrt(sum(tmp(h, 1:2).^2));
            for z_ref = -1:-1:-50
                p = make_arc(Tnodes_new, Tedges, Cz_new, tmp(h, :), [tmp(h, 1:2)*80/l z_ref]);
                p = p(1:end-1, :);
                if p(end, 3)<0
                    break
                end
            end
        end
        id = find(p(:, 3)<0, 1);
        X = p(id, :)+(p(id-1, :)-p(id, :))/(p(id, 3)-p(id-1, 3))*p(id, 3);
        try
            if X(2)>0
                p = make_arc(Tnodes_new, Tedges, Fpz_new, [Fpz_new(1) X(2) 0], X);
            else
                p = make_arc(Tnodes_new, Tedges, Fpz_new, Oz_new, X);
            end
            X_percent(h) = sum(sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2)))/X_total;
        catch % in case of concavity
            X_percent(h) = sqrt(sum((X-Fpz_new).^2, 2))/X_total;
        end
        p = make_arc(Tnodes_new, Tedges, Cz_new, [X(1:2) Cz_new(3)], X);
        Y_total = sum(sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2)));
        try
            if tmp(h, 3)<Cz_new(3) % reference below Cz
                p = make_arc(Tnodes_new, Tedges, Cz_new, [X(1:2) Cz_new(3)], tmp(h, :));
            else % reference above Cz
                p = make_arc(Tnodes_new, Tedges, Cz_new, [Cz_new(1:2) tmp(h, 3)], tmp(h, :));
            end
            Y_percent(h) = sum(sqrt(sum((p(1:end-1, :)-p(2:end, :)).^2, 2)))/Y_total;
        catch % when segment length is smaller than a triangle
            Y_percent(h) = sqrt(sum((X-Cz_new).^2, 2))/Y_total;
        end
    end
%     save('extended_brain_projection.mat', 'X_percent', 'Y_percent', '-append')
end
tet.X_percent = X_percent*100;
tet.Y_percent = Y_percent*100;
%% Write the table
% writetable(tet(:, [9 3:8 10:14]), 'tet_code2MNI_lookup_extended.xlsx', 'Sheet', 'Reference', 'WriteMode', 'overwritesheet')