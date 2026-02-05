addpath ./../src/ ; 
addpath ./../data/ ; 

%% Chargements de donnees pour CHOUCROUTE
load mini_wall_b.mat ; 
load mini_wall_b_gt.mat ; 
d_b = poissrnd(d) ; 

W = length(wavelengths) ; 
Ko = length(unique(lb_gt)) - 1 ;
s_ref_gt = zeros(W, Ko) ;
cube_gt_interp = reshape(cube_gt_interp, [], W) ; 
for k = 1 : Ko
    s_ref_gt(:, k) = mean(cube_gt_interp(lb_gt(:) == k, :)) ; 
end

%% CHOUCROUTE
tic
[im_lb, s_ref] = choucroute(d_b, H, 0.1, "order", "random", "save", "results_") ;
toc

% Affichage
[W, K] = size(s_ref) ; 
cmap = [0,0,0 ; rand(K, 3)] ; 
cmap_gt = [0,0,0 ; rand(Ko, 3)] ; 

figure ; 
% Carte des labels GT
subplot 221 ; imagesc(lb_gt) ; colormap(gca, cmap_gt) ; axis equal ; 
axis tight ; xticks([]) ; yticks([]) ; title("Labels GT") ; 
colorbar("Ticks", [0, Ko], "TickLabels", ["NC", num2str(Ko)]) ;
% Carte des labels CHOUCROUTE
subplot 222 ; imagesc(im_lb) ; colormap(gca, cmap) ; axis equal ; 
axis tight ; xticks([]) ; yticks([]) ; title("Labels CHOUCROUTE") ; 
colorbar("Ticks", [0, K], "TickLabels", ["NC", num2str(K)]) ; 
% Spectres de ref GT
subplot 223 ; hold on ; 
for k = 1 : Ko 
    plot(wavelengths, s_ref_gt(:, k), "Color", cmap_gt(k+1, :), "LineWidth", 2) ;
end ; hold off ; axis tight ; xlabel("Longueurs d'onde (nm)") ; 
ylabel("Amplitude") ; title("Spectres GT") ; ylim([0 2e3]) ; 
% SPectres de ref CHOUCROUTE
subplot 224 ; hold on ; 
for k = 1 : K 
    plot(wavelengths, s_ref(:, k), "Color", cmap(k+1, :), "LineWidth", 2) ;
end ; hold off ; axis tight ; xlabel("Longueurs d'onde (nm)") ; 
ylabel("Amplitude") ; title("Spectres estimes") ; ylim([0 2e3]) ; 
exportgraphics(gcf, "results.pdf", "ContentType", "vector") ; 

