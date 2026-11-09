clear;  close all;  clc;
% in limiti inserire guardando il wenglor da dove a dove si vedono i punti

% claibrazione 1: 
% limite_min = -9; limite_max = 26;
% cartella = fullfile('calibrazioni', 'CAL1/');

% calibrazione 2
% limite_min = -6; limite_max = 24;
% cartella = fullfile('calibrazioni', 'CAL2/');

%calibrazione 3
% limite_min = -6; limite_max = 21;
% cartella = fullfile('calibrazioni', 'CAL3/');

% calibrazione 4
%limite_min = -6; limite_max = 20;
%cartella = fullfile('calibrazioni', 'CAL4/');

% calibrazione 5
limite_min = -8; limite_max = 28;
cartella = fullfile('calibrazioni', 'CAL5/');

% calibrazione 6
% limite_min = -10; limite_max = 25;
% cartella = fullfile('calibrazioni', 'CAL6/');

% Assicurati che questi siano i numeri associati ai tuoi file!
indici_file = 0:1:19; 
n_misure = length(indici_file);

% Prealloco i vettori per salvare i centri di ogni singola posizione
centri_x = zeros(n_misure, 1);
centri_y = zeros(n_misure, 1);
diametri = zeros(n_misure, 1);

% --- ESTRAZIONE DEI CENTRI POSIZIONE PER POSIZIONE ---
for i = 1:n_misure
    k = indici_file(i);
    T = readmatrix(fullfile(cartella, sprintf('T%d.txt', k)));
    XY = T(:,1:2);
    
    % Filtro i punti
    XY = XY(~any(XY==0,2) & XY(:,1) >= limite_min & XY(:,1) <= limite_max, :);
    
    % Estraggo x e y SOLO per questa specifica acquisizione
    x = XY(:,1); 
    y = XY(:,2);
    
    % Fit del cerchio (Kåsa) applicato SOLO ai punti correnti
    p = [x y ones(size(x))]\(-(x.^2 + y.^2));    
    R = sqrt((p(1)^2 + p(2)^2) / 4 - p(3));

    % Salvo il centro trovato per questa acquisizione
    centri_x(i) = -p(1)/2;  
    centri_y(i) = -p(2)/2;
    diametri(i) = 2 * R;
end

% --- CALCOLO METROLOGICO DELL'ASSE E DELL'INCERTEZZA (Tipo A) ---

% L'asse di rotazione reale è la MEDIA della traiettoria dei centri
xc_asse = mean(centri_x);
yc_asse = mean(centri_y);
D_medio = mean(diametri);

% L'incertezza sull'asse è valutata statisticamente (Deviazione Standard)
u_xc = std(centri_x);
u_yc = std(centri_y);
% perchè non uso il fratto radice di n, ovvero
% u_xc = std(centri_x) / size(centri_x, 1);
% u_yc = std(centri_y) / size(centri_y, 1);
% Perché in questo caso non sto calcolando l'incertezza di ripetibilità su un misurando statico
% , ma sto quantificando l'ampiezza fisica dell'eccentricità del sistema di rotazione. L'ampiezza di questo sfarfallio è descritta dalla deviazione standard della popolazione dei centri
% , che è il parametro di dispersione corretto da inserire come perturbazione nel nostro modello Monte Carlo

% --- STAMPA DEI RISULTATI ---
fprintf('--- TARATURA ASSE ROTATORE ---\n');
fprintf('xc_asse = %.6f mm  (Incertezza u_xc = %.6f mm)\n', xc_asse, u_xc);
fprintf('yc_asse = %.6f mm  (Incertezza u_yc = %.6f mm)\n', yc_asse, u_yc);
fprintf('D_medio = %.6f mm\n', D_medio);

% Grafico per visualizzare l'eccentricità
figure('Color', 'w');
plot(centri_x, centri_y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Centri alle varie pose'); hold on;
plot(xc_asse, yc_asse, 'b+', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Asse Rotatore (Media)');
axis equal; grid on;
title('Taratura Asse Rotatore (Eccentricità Pezzo Bianco)');
xlabel('X Laser [mm]'); ylabel('Y Laser [mm]');
legend('Location','best');