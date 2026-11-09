clear; clc; close all;

%% 1) PARAMETRI DI CALIBRAZIONE (Incertezza del Setup)
% % calibrazione 1
% xc_nom = 16.050882; 
% yc_nom = 231.969615; 
% u_xc   = 0.096215;
% u_yc   = 0.103869;

% % calibrazione 2
% xc_nom = 22.224682; 
% yc_nom = 221.444876; 
% u_xc   = 0.077787;
% u_yc   = 0.158433;

% calibrazione 3
xc_nom = 24.999731; 
yc_nom = 212.893633; 
u_xc   = 0.077760;
u_yc   = 0.177220;

% calibrazione 4
% xc_nom = 26.097778; 
% yc_nom = 206.707625; 
% u_xc   = 0.026749;
% u_yc   = 0.119169;

R_sfera = 8.5 / 2; % Raggio della sfera per calcolo D0 (da CAD)

%% 2) SETUP MONTE CARLO
M = 30000; % 1000 per test veloce, 30000 solido 

% Generazione asse perturbato
xc_sim = xc_nom + randn(M, 1) * u_xc; % grazie a M genera vettore di M righe
yc_sim = yc_nom + randn(M, 1) * u_yc;

Diametri_In_MC  = zeros(M, 1);
Diametri_Out_MC = zeros(M, 1);
Diametri_D0_MC  = zeros(M, 1);

% errore meccanico del banco. Il pezzo sfarfalla. L'asse di rotazione varia 
% in modo casuale seguendo una distribuzione Gaussiana centrata in xc_nom, 
% con una dispersione (Incertezza di Tipo A) pari a u_xc
% Si generano 10.000 centri di rotazione fisicamente possibili

%% 3) PRE-CARICAMENTO DATI
cartella_rilevazione = fullfile('rilevazioni/rilevazioni_cal3/', '2026_04_23_RIL_2/');
angles = 0:30:330;
Dati_Grezzi = cell(length(angles), 1);

fprintf('Caricamento dati grezzi in corso...\n');
for i = 1:length(angles)
    v = angles(i);
    nome_file = fullfile(cartella_rilevazione, sprintf('%d.txt', v));
    if exist(nome_file, 'file')
        M_temp = readmatrix(nome_file);
        Dati_Grezzi{i} = M_temp(M_temp(:,1) ~= 0, 1:2);
    else
        warning('File %s non trovato.', nome_file);
    end
end
fprintf('Inizio simulazione Monte Carlo (%d iterazioni)...\n', M);

%% 4) MOTORE MONTE CARLO
for iter = 1:M
    xc = xc_sim(iter);
    yc = yc_sim(iter);
    
    % Allineamento
    distanza_sys_rif = hypot(xc, yc);
    phiX = atan2d(yc, xc);
    Ralign = [cosd(-phiX) -sind(-phiX); sind(-phiX) cosd(-phiX)];
    R_sys = distanza_sys_rif;

    P_temp_cell = cell(length(angles), 1); % cell array per i blocchi di punti

    for i = 1:length(angles)
        if isempty(Dati_Grezzi{i}), continue; end
        v = angles(i);

        % Trasformazioni geometriche
        P0 = (Ralign * Dati_Grezzi{i}(:,1:2).').';
        P  = ([cosd(v) -sind(v); sind(v) cosd(v)] * P0.').';
        O  = [R_sys * cosd(v), R_sys * sind(v)];

        % Salviamo i punti nel "cassetto" i-esimo della cell array
        P_temp_cell{i} = [P(:,1)-O(1), P(:,2)-O(2)];
    end

    % Concatena tutto in un colpo solo per maggiore velocità di esecuzione
    P_tot = vertcat(P_temp_cell{:});
    
    % --- 4.1 FIT KASA PER INTERNO/ESTERNO ---
    x = P_tot(:,1); y = P_tot(:,2);
    p_temp = [x y ones(size(x))]\(-(x.^2 + y.^2));   
    xc_t = -p_temp(1)/2; yc_t = -p_temp(2)/2;
    theta = atan2(y - yc_t, x - xc_t);      
    r     = hypot(x - xc_t, y - yc_t);

    nbins = 72; edges = linspace(-pi, pi, nbins+1);
    [b,~] = discretize(theta, edges);
    qLow = 0.10; qHigh = 0.90; w = 0.002;
    rootPts = []; tipPts = []; 

    for k_bin = 1:nbins
        idx = (b == k_bin);
        if nnz(idx) < 5, continue; end                 
        rk = r(idx); xk = x(idx); yk = y(idx);
        r_low = quantile(rk, qLow); r_high = quantile(rk, qHigh);
        rootPts = [rootPts; xk((rk>=r_low*(1-w)) & (rk<=r_low*(1+w))), yk((rk>=r_low*(1-w)) & (rk<=r_low*(1+w)))];
        tipPts  = [tipPts;  xk((rk>=r_high*(1-w)) & (rk<=r_high*(1+w))), yk((rk>=r_high*(1-w)) & (rk<=r_high*(1+w)))];
    end

    [~, ~, R_in]  = fit_circle_kasa(rootPts);
    [~, ~, R_out] = fit_circle_kasa(tipPts);
    Diametri_In_MC(iter)  = 2 * R_in;
    Diametri_Out_MC(iter) = 2 * R_out;

    % --- 4.2 FILTRAGGIO SELETTIVO DEI FIANCHI ---
    % Calcolo altezza del dente in questa iterazione
    tooth_depth = R_out - R_in;
    margin = tooth_depth * 0.15; % Taglio via il 15% dal fondo e 15% dalla cima (i raccordi)
    
    r_val = hypot(P_tot(:,1), P_tot(:,2));
    mask_keep = (r_val > R_in + margin) & (r_val < R_out - margin);
    P_keep = P_tot(mask_keep, :);

    % Se il rumore ha distrutto il profilo in questo giro e P_keep è vuoto, salta l'iterazione
    if size(P_keep,1) < 20
        Diametri_D0_MC(iter) = NaN;
        continue;
    end

    % --- 4.3 ESTRAZIONE RETTE (PCA) ---
    % è un fit delle rette attraverso la regressione lineare
    angoli_keep = atan2(P_keep(:,2), P_keep(:,1));
    [angoli_keep, sort_idx] = sort(angoli_keep);
    P_keep = P_keep(sort_idx, :);
    
    salti = diff(angoli_keep);
    salti(salti < 0) = salti(salti < 0) + 2*pi;
    idx_salti = find(salti > deg2rad(1.5));
    
    linee = []; 
    start_idx = 1;
    for g = 1:length(idx_salti)+1
        if g <= length(idx_salti), end_idx = idx_salti(g); else, end_idx = length(angoli_keep); end
        n_pts = end_idx - start_idx + 1;
        
        if n_pts >= 8 
            % Ultima pulizia: elimino i primissimi e ultimissimi punti del segmento
            trim = round(n_pts * 0.10); 
            pts = P_keep(start_idx+trim : end_idx-trim, :);
            if size(pts,1) >= 4
                centroide = mean(pts);
                coeff = pca(pts);
                linee = [linee; centroide, coeff(:,1)'];
            end
        end
        start_idx = end_idx + 1;
    end
    
    % --- 4.4 INCASTRO SFERE CORRETTO (CRESTA vs VALLE) ---
    Centri_Sfere = [];
    R_mid = (R_in + R_out) / 2; % SOGLIA DINAMICA INTELLIGENTE
    
    for k = 1:size(linee, 1)
        idx1 = k; idx2 = mod(k, size(linee,1)) + 1;
        C1 = linee(idx1, 1:2); v1 = linee(idx1, 3:4);
        C2 = linee(idx2, 1:2); v2 = linee(idx2, 3:4);
        
        A = [v1(1), -v2(1); v1(2), -v2(2)];
        B = [C2(1) - C1(1); C2(2) - C1(2)];
        if cond(A) > 1e10, continue; end 
        t = A \ B;
        V = C1 + t(1)*v1; 
        
        % CONTROLLO CRESTA vs VALLE DINAMICO
        if norm(V) > R_mid 
            continue; 
        end
        
        d1 = C1 - V; d1 = d1 / norm(d1);
        d2 = C2 - V; d2 = d2 / norm(d2);
        
        ub = d1 + d2; 
        if norm(ub) == 0, continue; end
        ub = ub / norm(ub);
        if dot(ub, V) < 0, ub = -ub; end 
        
        cos_2alpha = dot(d1, d2);
        sin_alpha = sqrt(max((1 - cos_2alpha)/2, 0.001));
        
        dist_Vs = R_sfera / sin_alpha;
        Cs = V + dist_Vs * ub;
        
        % Controllo di sanità geometrico
        if norm(Cs) > R_in && norm(Cs) < R_out + 2*R_sfera
            Centri_Sfere = [Centri_Sfere; Cs];
        end
    end
    
    if size(Centri_Sfere, 1) >= 3
        [~, ~, R_centers] = fit_circle_kasa(Centri_Sfere);
        Diametri_D0_MC(iter) = (2 * R_centers) + (2 * R_sfera);
    else
        R_centers = NaN; 
        Diametri_D0_MC(iter) = NaN;
    end
end

Diametri_D0_MC = Diametri_D0_MC(~isnan(Diametri_D0_MC));
M_valid = length(Diametri_D0_MC);

if M_valid == 0
    error('Nessuna iterazione valida trovata. Controlla il caricamento dati!');
end

%% 5) AGGIUNTA INCERTEZZA SENSORE LASER
Wenglor_Linearity_R  = 0.0325; 
Wenglor_Resolution_R = 0.014;  

dist_lin = 2 * (-Wenglor_Linearity_R + (2*Wenglor_Linearity_R) .* rand(M_valid, 1));
dist_res = 2 * (-(Wenglor_Resolution_R/2) + Wenglor_Resolution_R .* rand(M_valid, 1));

D_In_Fin  = Diametri_In_MC(1:M_valid)  + dist_lin + dist_res;
D_Out_Fin = Diametri_Out_MC(1:M_valid) + dist_lin + dist_res;
D_D0_Fin  = Diametri_D0_MC + dist_lin + dist_res;

% Pulizia outlier finali estremi dovuti alle aberrazioni PCA di Monte Carlo
mask_clean = (D_D0_Fin > 180 & D_D0_Fin < 192);
D_D0_Fin = D_D0_Fin(mask_clean);

%% 6) STATISTICA (Percentili 95%)
Mean_In  = mean(D_In_Fin);  U_In  = (prctile(D_In_Fin, 97.5) - prctile(D_In_Fin, 2.5)) / 2;
Mean_Out = mean(D_Out_Fin); U_Out = (prctile(D_Out_Fin, 97.5) - prctile(D_Out_Fin, 2.5)) / 2;
Mean_D0  = mean(D_D0_Fin);  U_D0  = (prctile(D_D0_Fin, 97.5) - prctile(D_D0_Fin, 2.5)) / 2;

%% 7) STAMPA IN COMMAND WINDOW
clc;
% Arrotonda le incertezze a 2 cifre significative (Regole GUM)
U_In_round  = round(U_In, 2, 'significant');
U_Out_round = round(U_Out, 2, 'significant');
U_D0_round  = round(U_D0, 2, 'significant');

fprintf('=== RISULTATI METROLOGICI DEFINITIVI (Iterazioni valide=%d) ===\n', length(D_D0_Fin));
% Stampa arrotondando in automatico anche i diametri allo stesso numero di decimali (.2f)
fprintf('D_in   = (%.2f +/- %.2f) mm (P=95%%)\n', Mean_In, U_In_round);
fprintf('D_out  = (%.2f +/- %.2f) mm (P=95%%)\n', Mean_Out, U_Out_round);
fprintf('D_0    = (%.2f +/- %.2f) mm (P=95%%)\n', Mean_D0, U_D0_round);
fprintf('===============================================================\n');

%% 8) PLOT dell'ultima Iterazione valida
figure('Name', 'Analisi Geometrica (Ultima Iterazione)', 'Color', 'w');

% 1. Profilo grezzo della ruota dentata (grigio)
plot(P_tot(:,1), P_tot(:,2), '.', 'Color', [0.8 0.8 0.8]); hold on; axis equal; grid on;

% 2. Fianchi segmentati (verde)
if ~isempty(P_keep)
    plot(P_keep(:,1), P_keep(:,2), '.g', 'MarkerSize', 8, 'DisplayName', 'Fianchi Segmentati');
end

% 3. Rette fittate con la PCA (magenta)
for k = 1:size(linee,1)
    C = linee(k, 1:2); v = linee(k, 3:4);
    p1 = C - v*8; p2 = C + v*8;
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'm-', 'LineWidth', 1.5, 'HandleVisibility','off');
end

% Aumento i punti a 100 per avere circonferenze disegnate in modo più fluido
theta = linspace(0, 2*pi, 100);

% 4. Sfere (arancioni per staccare dal grigio del profilo)
for k = 1:size(Centri_Sfere,1)
    plot(Centri_Sfere(k,1) + R_sfera*cos(theta), Centri_Sfere(k,2) + R_sfera*sin(theta), '-', 'Color', '#D95319', 'LineWidth', 2, 'HandleVisibility','off');
    plot(Centri_Sfere(k,1), Centri_Sfere(k,2), 'k+', 'MarkerSize', 5, 'HandleVisibility','off');
end

% 5. Diametri Interno ed Esterno
plot(R_in*cos(theta), R_in*sin(theta), 'b--', 'LineWidth', 1.5, 'DisplayName', 'D_{in}');
plot(R_out*cos(theta), R_out*sin(theta), 'r--', 'LineWidth', 1.5, 'DisplayName', 'D_{out}');

% 6. Traiettoria centri (sottile) e Diametro D_0 (tangente esterna, marcata)
if exist('R_centers', 'var') && ~isnan(R_centers)
    % Traiettoria dei centri delle sfere: alleggerita (sottile e tratteggiata)
    plot(R_centers*cos(theta), R_centers*sin(theta), 'k--', 'LineWidth', 0.5, 'DisplayName', 'Traiettoria Centri Sfere');
    
    % Calcolo e plot del D_0 (circonferenza tangente esternamente alle sfere)
    R_D0 = R_centers + R_sfera;
    plot(R_D0*cos(theta), R_D0*sin(theta), 'k-', 'LineWidth', 2.5, 'DisplayName', 'D_0 (Tangente Esterna Sfere)');
end

title('Analisi Geometrica - Ultima iterazione modello a sfere', 'FontSize', 14); 
legend('Location', 'best');
%% 9) ISTOGRAMMI CON LIMITI DI CONFIDENZA AL 95%

figure('Color', 'w', 'Name', 'Incertezze (P=95%) Monte Carlo Din'); 
histogram(D_In_Fin, 30, 'FaceColor', 'b'); hold on; grid on;
xline(Mean_In, '-k', 'LineWidth', 2, 'DisplayName', 'Media');
xline(prctile(D_In_Fin, 2.5), '--r', 'LineWidth', 2, 'DisplayName', 'LCL (2.5%)');
xline(prctile(D_In_Fin, 97.5), '--r', 'LineWidth', 2, 'DisplayName', 'UCL (97.5%)');
title('Incertezza D_{in}'); xlabel('Diametro [mm]'); ylabel('Frequenza'); legend('Location', 'best');
xticks(floor(min(D_In_Fin)) : 0.2 : ceil(max(D_In_Fin))); xtickformat('%.1f');

figure('Color', 'w', 'Name', 'Incertezze (P=95%) Monte Carlo Dout'); 
histogram(D_Out_Fin, 30, 'FaceColor', 'r'); hold on; grid on;
xline(Mean_Out, '-k', 'LineWidth', 2, 'DisplayName', 'Media');
xline(prctile(D_Out_Fin, 2.5), '--r', 'LineWidth', 2, 'DisplayName', 'LCL (2.5%)');
xline(prctile(D_Out_Fin, 97.5), '--r', 'LineWidth', 2, 'DisplayName', 'UCL (97.5%)');
title('Incertezza D_{out}'); xlabel('Diametro [mm]'); ylabel('Frequenza'); legend('Location', 'best');
xticks(floor(min(D_Out_Fin)) : 0.2 : ceil(max(D_Out_Fin))); xtickformat('%.1f');

figure('Color', 'w', 'Name', 'Incertezze (P=95%) Monte Carlo D esterno sfere'); 
histogram(D_D0_Fin, 30, 'FaceColor', 'g'); hold on; grid on;
xline(Mean_D0, '-k', 'LineWidth', 2, 'DisplayName', 'Media');
xline(prctile(D_D0_Fin, 2.5), '--r', 'LineWidth', 2, 'DisplayName', 'LCL (2.5%)');
xline(prctile(D_D0_Fin, 97.5), '--r', 'LineWidth', 2, 'DisplayName', 'UCL (97.5%)');
title('Incertezza D_{0} (Sfere)'); xlabel('Diametro [mm]'); ylabel('Frequenza'); legend('Location', 'best');
xticks(floor(min(D_D0_Fin)) : 0.2 : ceil(max(D_D0_Fin))); xtickformat('%.1f');

%% -----------------------------------------
%% FUNZIONI LOCALI MATEMATICHE
%% -----------------------------------------
function [xc, yc, R] = fit_circle_kasa(P)
    if size(P,1) < 3, xc=NaN; yc=NaN; R=NaN; return; end
    x = P(:,1); y = P(:,2);
    p = [x y ones(size(x))]\(-(x.^2 + y.^2));   
    xc = -p(1)/2; yc = -p(2)/2;
    R  = sqrt(max((p(1)^2 + p(2)^2)/4 - p(3), 0));
end