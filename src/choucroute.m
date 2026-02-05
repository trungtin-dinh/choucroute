function [im_lb, s_ref_ex] = choucroute(varargin)
% [im_lb, s_ref_ex] = CHOUCROUTE(d, H, T_psi, "PARAM1", PARAM1, "PARAM2", PARAM2, ...)
%
% Algorithme iteratif de classification non supervisee de donnees hyperspectrales 
% codees par tests statistiques en 3 etapes principales : detection,
% croissance et fusion.
%
% Entrees obligatoires : 
% d         (RxCxS)     : Donnees codees 
% H         (RxCxWxS)   : Cube de filtrage 
% T_psi     (1x1)       : Seuil de variation des coefficients de variabilite 
%                       spectrale intraclasse |1 - psi| 
% Sorties : 
% im_lb     (RxC)       : Image des labels de classification 
% s_ref     (WxK)       : Spectres de reference
%
% Entrees optionnelles : 
% pan       (RxC)       : Image panchromatique 
% T_pan     (1x1)       : Seuil d'intensite sur l'image panchromatique
%                       "auto" : choisi par methode d'Otsu
%                       Defaut : T_pan = 0
% P_d       (1x1)       : Taille du bloc a considerer pour la detection 
%                       (bloc de taille (2P+1) x (2P+1))  
%                        Defaut : P_d = ceil(ceil(sqrt(W/S))/2)
% P_c       (1x1)       : Taille du bloc a considerer pour la croissance 
%                       (bloc de taille (2P+1) x (2P+1))  
%                       Defaut : P_d = ceil(ceil(sqrt(W/S))/2)
% N_d       (1x1)       : Nombre de pixels minimum de la region homogene a
%                       l'etape de detection pour eviter qu'elle soit trop petite 
% N_f       (1x1)       : Nombre de pixels minimum de la region homogene a
%                       l'etape de fusion pour eviter qu'elle soit trop petite 
% masque    (RmxCm)     : Masque de croissance de region
% ordre     (string)    : Ordre de traitement de donnees 
%                       "none"      : spatial de gauche a droite, 
%                                   du haut en bas de l'image (par defaut)
%                       "ascend"    : croissant des intensites de pan
%                       "descend"   : decroissant des intensites de pan
%                       "random"    : aleatoire
%                       "custom"    : defini par l'utilisateur
% test      (string)    : Type de test 
%                       "sw"    : Shapiro-Wilk (par defaut)
%                       "ks"    : Kolmogorov-Smirnov
%                       "ad"    : Anderson-Darling
% alpha     (1x1)       : Niveau de signification des tests statistiques 
% est_psi_f (boolean)   : Estimation ou non des coefficients de variabilite 
%                       spectrale a l'etape de fusion des classes 
% mu        (1x1)       : Coefficient de regularisation pour l'estimation 
%                       quadratique du spectre 
% save      (boolean)   : Sauvegarder des sorties en un fichier .mat
% nom       (string)    : Nom du fichier de sauvegarde
%
% R         : Nombre de lignes
% C         : Nombre de colonnes
% W         : Nombre de longueurs d'onde
% S         : Nombre d'acquisition
% Rm, Cm    : Nombre de lignes et de colonnes du masque

    d = varargin{1} ; 
    H = varargin{2} ; 
    T_psi = varargin{3} ; 
    [R, C, W, S] = size(H) ;

    % Valeurs par defaut
    im_p = sum(d, 3) ; 
    T_pan = 0 ; 
    P_d = ceil(ceil(sqrt(W/S))/2) ;  
    P_c = P_d ; 
    N_d = ceil(W/S) ;
    N_f = 2*N_d ; 
    alpha = 0.05 ;     
    mask = [0 1 0 ; 1 1 1 ; 0 1 0] ; 
    ordre = "ascend" ; 
    test_type = "sw" ; 
    est_psi_f = false ; 
    sauvegarder = false ; 
    mu = [] ; 

    for n = 4:nargin
        if isstring(varargin{n})
            param = lower(varargin{n}) ; 
            switch param
                case "pan" ; im_p = varargin{n+1} ; 
                case {"pd", "p_d"} ; P_d = varargin{n+1} ; 
                case {"pc", "p_c"} ; P_c = varargin{n+1} ; 
                case {"nd", "n_d"} ; N_d = varargin{n+1} ; 
                case {"nf", "n_f"} ; N_f = varargin{n+1} ; 
                case {"mask", "masque"} ; mask = varargin{n+1} ; 
                case {"order", "ordre"}
                    ordre = lower(varargin{n+1}) ;  
                    if strcmp(ordre, "custom")
                        px_rest = varargin{n+2} ;                 
                    end                            
                case "test" ; test_type = varargin{n+1} ; 
                case "alpha" ; alpha = varargin{n+1} ;
                case "mu" ; mu = varargin{n+1} ;       
                case "est_psi_f" ; est_psi_f = varargin{n+1} ; 
                case "save" ; sauvegarder = true ; nom = varargin{n+1} ;
                case "t_pan"                    
                    if strcmpi("auto", varargin{n+1})
                        T_pan = graythresh(im_p(:)/max(im_p(:)))*max(im_p(:)) ;
                    else
                        T_pan = varargin{n+1} ; 
                    end            
            end
        end
    end

    % Reglage automatique du coefficient de regularisation SA
    if isempty(mu)
        mu = max(im_p(:))*2e-3 ; 
    end

    % Liste des pixels restants a traiter
    switch ordre
        case "ascend" ; [~, px_rest] = sort(reshape(im_p(im_p > T_pan), [], 1), ordre) ;
        case "descend" ; [~, px_rest] = sort(reshape(im_p(im_p > T_pan), [], 1), ordre) ;
        case "random" ; px_rest = 1:sum(im_p(:) > T_pan) ; px_rest = px_rest(randperm(length(px_rest))) ; 
        case "spatial" ; px_rest = 1:R*C ;
    end
    
    % Estimation du spectre de reference par NNLS
    f_H_temp = @(im_p, H_shape) reshape(permute(H_shape, [1, 3, 2]), [], W) ;
    f_s_ref = @(im_p, H_temp, D_ltD_l, d_h_ex) lsqnonneg(H_temp'*H_temp + mu*D_ltD_l, H_temp'*d_h_ex(:)) ;
    
    % Types de test
    switch test_type
        case "sw" ; f_test = @(d_e, d_p, alpha) swtest(d_e(:) - d_p(:), alpha) ;
        case "ks" ; f_test = @(d_e, d_p, alpha) kstest((d_e(:) - d_p(:))/std((d_e(:) - d_p(:))), alpha) ;
        case "ad" ; f_test = @(d_e, d_p, alpha) adtest((d_e(:) - d_p(:)), 'Alpha', alpha) ;
    end
     
    % Initialisation
    im_lb = zeros(R, C) ;       % Image des labels
    im_lb_old = ones(R, C) ;    % Image des labels de l'iteration precedente
    lb = 1 ;                    % Numero du label
    s_ref_ex = [] ; 

    % Redimensionner les variables pour faciliter l'extraction des donnees
    % dans les calculs qui suivent
    d_shape = reshape(d, [], S) ;       % Donnees codees
    H_shape = reshape(H, [], W, S) ;    % Cube de filtrage

    % Pre-calculer pour l'estimation des spectres
    % [1, -1, 0, ..., 0 ; ... ; 0, ..., 1, -1 ; 0, ..., 0, 1] (matrice de dimensions W x W)
    % [1, -1, 0, ..., 0 ; ... ; 0, ..., 1, -1] (matrice de dimensions W x (W-1))
    % Terme invariant dans les expressions de la boucle de calcul -> Calculer hors de la boucle pour reduire le cout de calcul
    D_l = diag(ones(W,1),0) + diag(-ones(W-1,1),1) ;
    D_l(end,:) = [] ;
    D_ltD_l = D_l'*D_l ;
     
    % Repeter l'algorithme jusqu'a tous les pixels soient traites, ou
    % l'image des labels n'evolue plus
    while (~isempty(px_rest)) && (sum(im_lb_old(:) - im_lb(:)) ~= 0) 
        im_lb_old = im_lb ;         % Resultats de l'iteration precedente 
        im_res = false(R, C) ;      % Positions des pixels a labeliser
    
        %%% Detection une region homogene 
        % Initialisation
        ind_h = 1 ;     % Indice du pixel central de la region homogene a tester dans la liste des pixels restants
        pval = 0 ;      % p-lavleur des tests
        while and((ind_h <= length(px_rest)), (pval < alpha))
            s_ref = zeros(W, 1) ; 
            pval = 0 ; 
            % Coordonnees spatiales qui correspondent a l'indice ind_h 
            [r_c, c_c] = ind2sub([R, C], px_rest(ind_h)) ;
            % Coordonnees des pixels autour de ce centre pour former un bloc
            r_h = max(1, r_c-P_d) : min(R, r_c+P_d) ;
            c_h = max(1, c_c-P_d) : min(C, c_c+P_d) ;
            [r_grid, c_grid] = meshgrid(r_h, c_h) ; 
            % Comparer ces coordonnees avec la liste des pixels restants pour eviter de doublement
            % px_h : indices des pixels de la region de depart a tester 
            px_h = intersect(sub2ind([R, C], r_grid(:), c_grid(:)), px_rest) ;         
            px_h(abs(1-im_p(px_h)/mean(im_p(px_h))) > T_psi) = [] ;
            % Si le nombre de pixels de la region a tester reste suffisant, 
            % alors estimer le spectre de reference de la region            
            if (length(px_h) > N_d) 
                % Estimation du spectre de reference par NNLS                 
                s_ref = f_s_ref(im_p(px_h), f_H_temp(im_p(px_h), H_shape(px_h, :, :)), D_ltD_l, d_shape(px_h, :)) ;
                s_ref(isnan(s_ref)) = 0 ;
                if (sum(s_ref) ~= 0)
                    % Carte des psi
                    im_psi = im_p/sum(s_ref) ;      
                    % Donnees codees extraites
                    d_h_o = d_shape(px_h, :) ;
                    % Prediction des donnees codees                        
                    d_h_p = squeeze(sum(reshape(im_psi(px_h).*s_ref', [length(px_h), W]).*H_shape(px_h, :, :), 2)) ; 
                    % Test statistiques
                    if ~ttest(d_h_o(:) - d_h_p(:), [], 'Alpha', alpha)
                        [~, pval] = f_test(d_h_o, d_h_p, alpha) ;
                    end
                end                    
            end
            ind_h = ind_h + 1 ; 
        end
        if (sum(s_ref) == 0), break ; end
        im_res(px_h) = true ;         
        
        % Supprimer les indices des pixels traites 
        px_rest = setdiff(px_rest, px_h, 'stable') ;  
        
        %%%%%%%%% Croissance de region
        % Tant qu'il reste encore de changements, continuer de croitre la 
        % region homogene actuelle
        im_res_old = false(R, C) ;
        while (sum(double(im_res_old(:)) - double(im_res(:))) ~= 0)
            im_res_old = im_res ;
            % ind_d : indices des pixels voisins de la region homogene 
            ind_d = find(imdilate(im_res, mask) - im_res) ;
            % Negliger les pixels ayant |1 - psi| > seuil_psi
            ind_d(abs(im_psi(ind_d)-1) > T_psi) = [] ;       
            if isempty(ind_d) ; break ; end
            % Pour chacun de ces pixels, prendre un bloc de taille 
            % (2P+1) autour et tester
            % Convertion de ces indices en coordonnees spatiales                    
            [r_c, c_c] = ind2sub([R, C], ind_d) ;
                        
            % Coordonnees des extremites des blocs correspondants a 
            % chaque pixel        
            r_first = max(1, r_c-P_c) ; r_end = min(R, r_c+P_c) ; 
            c_first = max(1, c_c-P_c) ; c_end = min(C, c_c+P_c) ; 
            % Attention, la taille des blocs peuvent etre 
            % differentes donc on ne peut pas les traiter en meme temps
            % On test bloc par bloc
            for ind_t = 1 : length(ind_d)            
                mat_ind = false(R, C) ;                 
                [r_grid, c_grid] = meshgrid(r_first(ind_t):r_end(ind_t), c_first(ind_t):c_end(ind_t)) ; 
                mat_ind(sub2ind([R, C], r_grid(:), c_grid(:))) = true ;
                % Negliger les pixels ayant |1 - psi| > seuil_psi, et les
                % pixels de la region homogene car sinon 2 fois les memes 
                % valeurs des residus pour la partie commune (histogramme)
                mat_ind(abs(1 - im_psi) > T_psi | im_res) = false ;
                mat_ind(~px_rest) = false ;
                mat_ind(mat_ind & (im_lb > 0)) =false ; 
                
                % S'il ne reste aucun pixel a tester, on passe a
                % l'iteration suivante
                if nnz(mat_ind)
                    % Extraction des donnees
                    d_o = d_shape(mat_ind, :) ;        % Donnees codees
                    H_ex = H_shape(mat_ind, :, :) ;     % Cube de filtrage                    
                    % Prediction des donnees codees
                    o_pre = im_psi(mat_ind).*s_ref' ; 
                    d_pre = squeeze(sum(o_pre.*H_ex, 2)) ; 
                    % Concatenation des donnees : 
                    % Si la region homogene de depart est plus grande que 
                    % le bloc a tester, on choisit aleatoireement quelques 
                    % pixels de la region homogene
                    if (size(d_h_o, 1) > size(d_o, 1))
                        idx = randperm(size(d_h_o, 1), size(d_o, 1)) ; 
                        d_ext_c = [reshape(d_h_o(idx, :), [], 1) ; d_o(:)] ;
                        d_pre_c = [reshape(d_h_p(idx, :), [], 1) ; d_pre(:)] ;
                    elseif (size(d_h_o, 1) < size(d_o, 1))
                        idx = randperm(size(d_o, 1), size(d_h_o, 1)) ; 
                        d_ext_c = [reshape(d_o(idx, :), [], 1) ; d_h_o(:)] ;
                        d_pre_c = [reshape(d_pre(idx, :), [], 1) ; d_h_p(:)] ;
                    else
                        d_ext_c = [d_h_o(:) ; d_o(:)] ;
                        d_pre_c = [d_h_p(:) ; d_pre(:)] ;
                    end
        
                    % Test statistiques
                    if ~ttest(d_ext_c(:) - d_pre_c(:), [], 'Alpha', alpha)
                        [~, pval] = f_test(d_ext_c, d_pre_c, alpha) ;
                    end
                    
                    % Labelisation des nouveaux pixels -> region elargie 
                    if (pval > alpha)
                        mat_ind(~px_rest) = false ;
                        im_res(mat_ind) = true ;                        
                        % px_rest = setdiff(px_rest, find(im_res), 'stable') ; 
                    end  
                end 
            end
        end
                        
        im_lb(im_res) = lb ;
        % Supprimer les pixels traites de la liste ind_rest
        px_rest = setdiff(px_rest, find(im_res), 'stable') ; 
        
        %%%%%%%%% Fusion des labels 
        % Concatenation des donnees de la region detectee et une des
        % regions precedemment detectee deux a deux
        new_label = true ; 
        for ind_lb = 1 : lb-1
            d_o = d_shape(im_res, :) ;
            d_lb = d_shape(im_lb == ind_lb, :) ; 
            
            if (length(d_o(:)) > length(d_lb(:)))
                idx = find(im_res) ; 
                idx = idx(randperm(length(idx), size(d_lb, 1))) ; 
                d_o = d_shape(im_res(idx), :) ; 
                d_lb = d_shape(im_lb == ind_lb, :) ;                 
                mat_ind = (im_lb == ind_lb) ; 
                mat_ind(idx) = true ; 
            elseif (length(d_o(:)) < length(d_lb(:)))           
                idx = find(im_lb(:) == ind_lb) ; 
                idx = idx(randperm(length(idx), size(d_o, 1))) ; 
                temp = im_lb == ind_lb ; 
                d_lb = d_lb(temp(idx), :) ;   
                mat_ind = im_res ; 
                mat_ind(idx) = true ; 
            else
                mat_ind = im_res ; 
                mat_ind(im_lb == ind_lb) = true ; 
            end
            d_c_o = [d_o ; d_lb] ; 
            
            % Comme a l'etape de detection, verifier d'abord si la plupart des pixels
            % sont suffisamment homogenes en intensite (au moins 2N pixels
            % car 2 classes), alors les calculs se poursuivent
            if (nnz(abs(1-im_p(mat_ind)/mean(im_p(mat_ind), "all")) < T_psi) > N_f)
                % Estimation du spectre de reference non normalise                        
                s_ref_c = f_s_ref(im_p(mat_ind), f_H_temp(im_p(mat_ind), H_shape(mat_ind, :, :)), D_ltD_l, d_c_o) ;            
                s_ref_c(isnan(s_ref_c)) = 0 ;
                % Prediction des donnees codees
                if est_psi_f 
                    psi_c = im_p(mat_ind)/sum(s_ref_c) ;
                else
                    psi_c = ones(size(im_p(mat_ind))) ;
                end
                d_c_p = squeeze(sum((psi_c.*s_ref_c').*H_shape(mat_ind, :, :), 2)) ; 
                % Tests sur les residus entre les donnees reelles et les 
                % donnees predites
                if ~ttest(d_c_o(:) - d_c_p(:), [], 'Alpha', alpha)
                    [~, pval] = f_test(d_c_o, d_c_p, alpha) ;
                end
                
                % Si le test est valide, alors on attribue le label
                % correspondant a la nouvelle region detectee
                if (pval > alpha)
                    im_lb(im_res) = ind_lb ;
                    new_label = false ;
                    s_ref_ex(:, ind_lb) = s_ref_c(:) ; 
                    break ;
                end
            end 
        end
        
        % Si une nouvelle region est detectee 
        % -> Incrementer le numero du label suivant
        if new_label
            im_lb(im_res) = lb ; 
            lb = lb+1 ; 
            s_ref_ex = [s_ref_ex, s_ref(:)] ; 
        end
    end         
        
    if sauvegarder
        full_nom = "Tpsi_" + num2str(T_psi) ; 
        for n = 4:nargin
            if isstring(varargin{n})
                param = lower(varargin{n}) ; 
                switch param                
                    case {"pd", "p_d"} ; full_nom = full_nom + "_" + param + "_" + num2str(P_d) ; 
                    case {"pc", "p_c"} ; full_nom = full_nom + "_" + param + "_" + num2str(P_c) ; 
                    case {"nd", "n_d"} ; full_nom = full_nom + "_" + param + "_" + num2str(N_d) ; 
                    case {"nf", "n_f"} ; full_nom = full_nom + "_" + param + "_" + num2str(N_f) ; 
                    case {"mask", "masque"} ; full_nom = full_nom + "_" + param + "_" + mask ; 
                    case {"order", "ordre"} ; full_nom = full_nom + "_" + param + "_" + ordre ;                     
                    case "test" ; full_nom = full_nom + "_" + param + "_" + test ; 
                    case "alpha" ; full_nom = full_nom + "_" + param + "_" + num2str(alpha) ; 
                    case "mu" ; full_nom = full_nom + "_" + param + "_" + num2str(mu) ;
                    case "est_psi_f" ; full_nom = full_nom + "_" + param + "_" + num2str(est_psi_f) ;                 
                    case "t_pan" ; full_nom = full_nom + "_T_" + num2str(T_pan) ;
                end
            end
        end
        s_ref = s_ref_ex ; 
        K = length(unique(im_lb)) - 1 ; 
        s_ref = zeros(W, K) ; 
        for k = 1 : K
            s_ref(:, k) = f_s_ref(im_p(im_lb == k), f_H_temp(im_p(im_lb == k), H_shape(im_lb(:) == k, :, :)), D_ltD_l, d_shape(im_lb(:) == k, :)) ;  
        end
        save(nom+full_nom+".mat", "im_lb", "s_ref") ; 
    end    
end

function [T, p, z] = swtestlite(X, alpha)
% SW_LITE Test de Shapiro-Wilk 
% Ce test de Shapiro-Wilk teste l'hypothese nulle (H0) selon laquelle un 
% echantillon est issu d'une population normalement distribuee
%
%   [T, p] = SW_LITE(X, alpha)
%   SW_LITE(X)
%
%   Entrees : 
%       X       :   Matrice de P variables a tester, chaque variable 
%                   contient N observations (N x P)
%                   Attention : N > 2
%       alpha   :   Niveau de signification (1 x 1) (optionnel)
%   Sorties : 
%       T       :   Ensemble des resultats des tests sur P variables (P x 1) 
%                   T = 0 : On ne rejete pas l'hypothese nulle (test valide)
%                   T = 1 : On rejete l'hypothese nulle
%                   H0 n'est pas rejetee quand p-value > alpha.
%       p       :   p-value (P x 1)
%       z       :   z-score (P x 1)
%
%   Valeur par defaut : 
%       alpha = 0.05

    % Statistique theorique de test pour 1 variable x : 
    %           sum(a*x)^2
    %   W = ------------------
    %       sum(x - mean(x))^2
    %
    % En pratique, differents calculs seront effectues en fonction du
    % nombre d'observations N (N > 2). On a generalement deux cas :   
    %   2 < N < 12 : La statistique du test se calcule independamment aux
    %   valeurs des observations.
    %   N >= 12 : La statistique du test se calcule en fonction des
    %   observations.    
    %
    % Les coefficients a sont calcules en fonction de m - les esperances 
    % des statistiques d'ordre d'un echantillon de variables iid suivant 
    % une loi normale.

    if (nargin < 2)
        alpha = 0.05 ;     
    end
    
    [N, P] = size(X) ; 
    X(isnan(X)) = 0 ;      
    
    m = norminv(((1:N)'-0.375)/(N+0.25));
    ms = sumsqr(m);

    % Pour la plupart des valeurs de N, les deux derniers coefficients sont
    % souvent calcules independamment.    
    u = 1/sqrt(N) ; 
    a_N        = (-2.706056*u^5) + (4.434685*u^4) - (2.07119*u^3) - (0.147981*u^2) + (0.221157*u) + (m(end,1)/sqrt(ms)) ;
    a_N_1  = (-3.582633*u^5) + (5.682633*u^4) - (1.752461*u^3) - (0.293762*u^2) + (0.042981*u) + (m(end-1,1)*(ms^-0.5)) ;
    
    % Calculs des coefficients a
    Phi = (ms - (2*m(end,1)^2) - (2*m(end-1,1)^2)) / (1 - (2*a_N^2) - (2*a_N_1^2)) ;
    if((N > 3) && (N < 6))
        Phi = (ms - 2*m(end,1)^2) / (1 - 2*a_N^2) ;
    end
    a = m/sqrt(Phi) ;  
    
    % Affectations des valeurs aux deux derniers coefficients a
    a(end,1) = a_N ;            % Coefficient N
    if(N >= 6)
        a(end-1,1) = a_N_1 ;    % Coefficient N-1
    end
    a(1,1) = -a(end,1) ;
    a(2,1) = -a(end-1,1) ;
    
    % Calculs des statistiques (moyenne et ecart-type) de la statistique 
    % du test Shapiro-Wilk (W)
    if(N == 3)
        a(3,1) = .70710678 ;
        a(1,1) = -a(end,1) ;
        a(2,1) = 0 ;
        W_mean = 0 ; 
        W_std = 1 ;
    elseif((N > 3) && (N < 12))        
        W_mean = (-0.0006714*N^3) + (0.025054*N^2) - (0.39978*N) + 0.544 ;
        W_std = exp( (-0.0020322*N^3) + (0.062767*N^2) - (0.77857*N) + 1.3822 ) ;  
        gamma1 = (0.459*N) - 2.273 ;	
    elseif(N >= 12)
        logN = log(N) ;         
        W_mean = (0.0038915*(logN)^3) - (0.083751*(logN)^2) - (0.31082*(logN)) - 1.5861 ;
        W_std = exp( (0.0030302*(logN^2)) - (0.082676*(logN)) - 0.4803 ) ;    
    end    
    
    % Calculs de p-value pour chaque variable de N observations 

    x = sort(X) ; 
    % Si on considere [x, a] comme une matrice (N x (P+1)), corrcoef
    % calcule les coefficients de correlation entre (P+1) colonnes. Les P
    % premieres lignes de la derniere colonne correspond aux coefficients 
    % de correlation entre chaque colonne de x et a.
    w_corr = corrcoef([x, a]) ;         % Coefficients de correlation
    W = w_corr(1:P, end).^2 ;           % Statistique du test
                    
    if(N == 3)                        
        p = 1.909859321 * (asin(sqrt(W))- 1.047198) ;
        z = norminv(p, W_mean, W_std) ; 
    elseif(N > 3) && (N < 12)            
        gamma2 = -log(1 - W) ;
        gamma3 = -log(gamma1 + gamma2) ;
        z = (gamma3-W_mean)/W_std ;   	
        p = 1 - normcdf(z) ;            % p-value
    else
        gW = log(1-W) ;                 % The transformation(g) value
        z = (gW-W_mean)/W_std ;         
        p = 1 - normcdf(z) ;            % p-value
    end    
    
    T = (p >= alpha) ;                  % On ne rejete pas l'hypothese nulle si alpha p-value > alpha
end
