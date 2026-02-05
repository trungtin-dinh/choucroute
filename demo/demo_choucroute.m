% Minimal demonstration of the CHOUCROUTE algorithm on synthetic data
% This script is self-contained and does not require external datasets.
addpath ./../src/ ; 

clear; clc; close all;

%% Dimensions (small on purpose: fast execution)
R = 30 ;    % Number of rows
C = 30 ;    % Number of columns
W = 50 ;    % Number of spectral bands
S = 3 ;     % Number of coded acquisitions
rng(0) ;    % For reproducibility

% Synthetic label map
im_gt = zeros(R, C) ;
im_gt(:, 1:10) = 1 ; 
im_gt(:, 11:20) = 2 ; 
im_gt(:, 21:30) = 3 ; 

% Synthetic reference spectra (ground truth)
s_ref_gt = zeros(W, 3) ; 
s_ref_gt(:, 1) = abs(randn(W, 1)) ;
s_ref_gt(:, 2) = abs(randn(W, 1)) + 50 ;
s_ref_gt(:, 3) = abs(randn(W, 1)) + 100 ;

% Synthetic filtering cube
H = cubefiltrage(W, dmd_fct(R, C, W, S)) ; 

% Intra-class variability coefficients
psi = 1 + 0.1*randn(R, C);

% Generate hyperspectral cube
cube_HS = zeros(R, C, W) ; 
for w = 1 : W
    cube_HS(:, :, w) = psi.*reshape(s_ref_gt(w, im_gt(:)), R, C) ;
end

% Generate coded measurements
d = poissrnd(squeeze(sum(cube_HS.*H, 3))) ; 

%% Run CHOUCROUTE
T_psi = 0.1 ;
[im_lb, s_ref] = choucroute(d, H, T_psi) ;

% Display results
figure ;
subplot 221 ; imagesc(im_gt) ; axis equal tight ; colorbar ;
title("Ground-truth labels") ;

subplot 222 ; imagesc(im_lb) ; axis equal tight ; colorbar;
title('CHOUCROUTE labels');

subplot 223 ; plot(s_ref_gt, "LineWidth", 2) ; axis tight ; 
title("Estimated reference spectra") ;

subplot 224 ; plot(s_ref, "LineWidth", 2) ; axis tight ; 
title("Estimated reference spectra") ;

%% Functions for filtering cube generating
function DMD = dmd_fct(R, C, W, S)        
    K = C + W - 1 ;                             % Number of columns of H
    H_model = zeros(R, W, S) ;                  % Initialize the DMD matrix
    
    for r = 1:R
        for n = 1:ceil(W/S)                     % Each portion of S columns
            dispo = (n-1)*S+1 : min(n*S, W) ;   % Ensure dispo doesn't exceed W
            for s = 1:min(S, length(dispo))     % Each acquisition (avoid exceeding dispo length)
                indice = randi([1, length(dispo)], 1) ; % Random index
                H_model(r, dispo(indice), s) = 1 ;      % Activate mirror
                dispo(indice) = [] ;                    % Remove used position
            end
        end
        
        % Handle remaining wavelengths
        dispo = (ceil(W/S) * S) + 1 : W ;               % Correct range
        for s = 1:min(S, length(dispo))                 % Avoid exceeding available positions
            indice = randi([1, length(dispo)], 1) ;
            H_model(r, dispo(indice), s) = 1 ;
            dispo(indice) = [] ;
        end
    end
    
    % Replicate W-periodically to reach K columns
    num_repeats = ceil(K/W) ;
    DMD = repmat(H_model, [1, num_repeats, 1]) ;
    
    % Truncate extra columns beyond K
    DMD = DMD(:, 1:K, :) ;
end

function Hcube = cubefiltrage(W, DMD, param_optique)
% CUBEFILTRAGE Generer le cube de filtrage
% ------------------------------
%	Hcube = CUBEFILTRAGE(W, DMD, param_optique)
%	Hcube = CUBEFILTRAGE(W, DMD)
% ------------------------------
%       Entrees :     
%           W : Nombre de longueurs d'onde (1)
%       	DMD : Matrice de configuration du DMD (R x (C+W-1) x N)
%         	param_optique : Paramtres de l'instrument (6)
%           	w0 : Longueur d'onde centrale (1)
%               x0_DMD : Centre du DMD (1)
%            	x0_objet : Centre de l'onjet sur le plan DMD (1)
%              	delta : Largeur de fente du DMD (1)
%              	alpha : Coefficient de dispersion (1)
%               beta : Facteur de grossissement (1)
%    	Sortie :      
%           Hcube : Cube de filtrage (R x C x W x N) 

    % Valeurs par defaut
    delta = 1 ;                                 % Largeur de fente du DMD 
    alpha = 1 ;                                 % Indice de dispersion 
    beta = 1 ;                                  % Facteur de grossissement 
    
    [R, K, N] = size(DMD) ;                     % Dimensions de la matrice du DMD
    C = K + 1 - W ;                             % Nombre de colonnes
    Hcube = zeros(R, C, W, N) ;                 % Initialisation
    
    if (nargin == 2)                            % Hcube = cubefiltrage(W, DMD)
        param_optique = [(W+1)/2, (W+C-1)/2, (W+C-1)/2, delta, alpha, beta] ; 
    end

    decalage = param_optique(2) + ((1:W) - param_optique(1))*param_optique(5)*param_optique(4) ;	% Decalage optique 
    debut = floor(decalage - C/2) + 1 ;         % Debut de la partie de la matrice du DMD a mettre dans le cube de filtrage
    fin = rem(debut + C - 1, K) ;               % Fin de la partie de la matrice du DMD a mettre dans le cube de filtrage obtenue par le reste de la division entre la position C apres le debut et K
    fin(fin == 0) = K ;                         % Si le reste de la division est nulle, la fin est egale a K

    for n = 1 : N                               % A chaque aquisition, on genere un cube de filtrage
        for w = 1 : W                           % A chaque longueur d'onde, on recupere une partie de la matrice du DMD et la mettre dans la cube de filtrage
            if (fin(w) - debut(w) + 1 < C)      % Si la distance entre le debut et la fin est plus petite que le nombre de colonnes du cube de filtrage
                Hcube(:, 1:(fin(w) - debut(w) + 1), w, n) = squeeze(DMD(:, debut(w) : K, n)) ;
                Hcube(:, fin(w) + 1 : end, w, n) = squeeze(DMD(:, 1 : fin(w), n)) ;
            else                                % Si la distance entre le debut et la fin est egale au nombre de colonnes du cube de filtrage          
                Hcube(:, :, w, n) = squeeze(DMD(:, debut(w) : fin(w), n)) ;
            end
        end
    end
end

