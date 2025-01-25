% ANNUALIZED
clear 
close all
clc

% Caricamento dei dati da Excel
data = readtable('Portfolio with youtube video.xlsx', 'Sheet', 'Sheet4');

% Estrarre le etichette degli asset
LABELS = data.Properties.VariableNames(2:11);

% Estrarre le date e i ritorni
DATES = data{:, 1}; % La prima colonna sono le date
RETURNS = data{:, 2:11}; % Le colonne restanti sono i ritorni degli asset

ASSET = size(RETURNS, 2); % Numero di asset
WINDOW_SIZE =  3 * 12; % Finestra temporale di 3 anni (in mesi)
NUM_WINDOWS = size(RETURNS, 1) - WINDOW_SIZE + 1; % Numero di finestre temporali

% Supponiamo di avere un tasso privo di rischio definito
Rf = 0.01/12; % esempio: 1%

% Inizializza matrici per i risultati
ROLLING_WEIGHTS = zeros(NUM_WINDOWS, ASSET);
ROLLING_RETURNS = zeros(NUM_WINDOWS, 1);
ROLLING_RISKS = zeros(NUM_WINDOWS, 1);

for w = 1:NUM_WINDOWS
    % Estrai i dati per la finestra temporale corrente
    WINDOW_RETURNS = RETURNS(w:w + WINDOW_SIZE - 1, :);
    
    % Calcola i ritorni attesi e la matrice di covarianza
    mean_returns = mean(WINDOW_RETURNS); % Ritorni medi mensili
    COV = cov(WINDOW_RETURNS); % Covarianza mensile
    
    % Inizializza variabili per il miglior portafoglio
    max_sharpe = -inf; 
    best_weights = zeros(1, ASSET);
    best_return = 0;
    best_risk = 0;
    
    % Esegui il loop su pesi casuali o generati in modo sistematico
    num_portfolios = 10000; % numero di portafogli da simulare
    for i = 1:num_portfolios
        % Genera pesi casuali per gli asset
        WTS = rand(1, ASSET);
        WTS = WTS / sum(WTS); % Normalizza i pesi per sommare a 1
        
        % Calcola rendimento e rischio per il portafoglio corrente
        port_return = sum(WTS .* mean_returns); % rendimento atteso mensile
        port_risk = sqrt(WTS * COV * WTS') * sqrt(12); % rischio annualizzato

        % Calcola lo Sharpe Ratio
        sharpe_ratio = (port_return - Rf) / port_risk;

        % Verifica se questo è il miglior portafoglio
        if sharpe_ratio > max_sharpe
            max_sharpe = sharpe_ratio;
            best_weights = WTS;
            best_return = port_return; % Ritorno mensile
            best_risk = port_risk; % Rischio annualizzato
        end
    end
    
    % Memorizza i risultati per il miglior portafoglio nella finestra corrente
    ROLLING_WEIGHTS(w, :) = best_weights;
    ROLLING_RETURNS(w) = best_return * 12; % Ritorno annualizzato
    ROLLING_RISKS(w) = best_risk;
end

% Grafico per la Prima Analisi (rendimenti e rischi variabili)
figure;
subplot(3,1,1);
plot(DATES(WINDOW_SIZE:end), ROLLING_RETURNS * 100);
title('Expected return of the portfolio (annualized)');
xlabel('Year');
ylabel('Expected Return (%)');

subplot(3,1,2);
plot(DATES(WINDOW_SIZE:end), ROLLING_RISKS * 100);
title('Risk of the portfolio (annualized)');
xlabel('Year');
ylabel('Risk (%)');

subplot(3,1,3);
plot(DATES(WINDOW_SIZE:end), sum(ROLLING_WEIGHTS(:, 1:7), 2) * 100, 'DisplayName', 'Equity');
hold on;
plot(DATES(WINDOW_SIZE:end), sum(ROLLING_WEIGHTS(:, 8:10), 2) * 100, 'DisplayName', 'Bond');
hold off;
title('Portfolio composition in time (annualized)');
xlabel('Year');
ylabel('Weight (%)');
legend('show');

% Focus sugli ultimi 5 anni (dal 2018 al 2023) - Prima Analisi
start_date = datetime(2018,1,1);
end_date = datetime(2023,12,31);
focus_indices = find(DATES(WINDOW_SIZE:end) >= start_date & DATES(WINDOW_SIZE:end) <= end_date);

if ~isempty(focus_indices)
    FOCUS_RETURNS = ROLLING_RETURNS(focus_indices);
    FOCUS_RISKS = ROLLING_RISKS(focus_indices);
    FOCUS_WEIGHTS = ROLLING_WEIGHTS(focus_indices, :);
    FOCUS_DATES = DATES(WINDOW_SIZE:end);
    FOCUS_DATES = FOCUS_DATES(focus_indices);

    figure;
    subplot(3,1,1);
    plot(FOCUS_DATES, FOCUS_RETURNS * 100);
    title('Expected return (2018-2023, annualized)');
    xlabel('Year');
    ylabel('Expected Return (%)');

    subplot(3,1,2);
    plot(FOCUS_DATES, FOCUS_RISKS * 100);
    title('Risk (2018-2023, annualized)');
    xlabel('Year');
    ylabel('Risk (%)');

    subplot(3,1,3);
    plot(FOCUS_DATES, sum(FOCUS_WEIGHTS(:, 1:7), 2) * 100, 'DisplayName', 'Equity');
    hold on;
    plot(FOCUS_DATES, sum(FOCUS_WEIGHTS(:, 8:10), 2) * 100, 'DisplayName', 'Bond');
    hold off;
    title('Portfolio composition (2018-2023, annualized)');
    xlabel('Year');
    ylabel('Weight (%)');
    ylim([0 100]); % Limita l'asse Y a 100
    legend('show');
end

% Seconda Analisi Modificata: Rendimento costante, rischio variabile (in base alla correlazione)
mean_exp_ret = mean(mean(RETURNS)); % Rendimento medio mensile

ROLLING_WEIGHTS_CONST_RETURN = zeros(NUM_WINDOWS, ASSET);
ROLLING_RISKS_CONST_RETURN = zeros(NUM_WINDOWS, 1); % Rischio variabile
ROLLING_RETURNS_CONST_RETURN = repmat(mean_exp_ret * 12, NUM_WINDOWS, 1); % Rendimento annualizzato

for w = 1:NUM_WINDOWS
    % Estrai i dati per la finestra temporale corrente
    WINDOW_RETURNS = RETURNS(w:w + WINDOW_SIZE - 1, :);
    
    % Matrice di covarianza annualizzata
    COV = cov(WINDOW_RETURNS); % Covarianza mensile

    % Inizializza variabili per il miglior portafoglio
    max_sharpe = -inf; 
    best_weights = zeros(1, ASSET);
    best_risk = 0;
    
    % Esegui il loop su pesi casuali o generati in modo sistematico
    num_portfolios = 10000; 
    for i = 1:num_portfolios
        % Genera pesi casuali per gli asset
        WTS = rand(1, ASSET);
        WTS = WTS / sum(WTS); 
        
        % Mantieni il rendimento costante
        port_return = mean_exp_ret * 12; % Rendimento annualizzato

        % Calcola il rischio per il portafoglio corrente (varia in base alla covarianza)
        port_risk = sqrt(WTS * COV * WTS') * sqrt(12); % Rischio annualizzato

        % Calcola lo Sharpe Ratio
        sharpe_ratio = (port_return - Rf * 12) / port_risk; % Rf annualizzato

        % Verifica se questo è il miglior portafoglio
        if sharpe_ratio > max_sharpe
            max_sharpe = sharpe_ratio;
            best_weights = WTS;
            best_risk = port_risk;
        end
    end
    
    % Memorizza i risultati per il miglior portafoglio nella finestra corrente
    ROLLING_WEIGHTS_CONST_RETURN(w, :) = best_weights;
    ROLLING_RISKS_CONST_RETURN(w) = best_risk;
end

% Grafico per la Seconda Analisi Modificata (rendimenti costanti, rischi variabili)
figure;
subplot(2,1,1);
plot(DATES(WINDOW_SIZE:end), ROLLING_RISKS_CONST_RETURN * 100);
title('Risk of the portfolio (constant return, variable risk, variable correlation, annualized)');
xlabel('Year');
ylabel('Risk (%)');

subplot(2,1,2);
plot(DATES(WINDOW_SIZE:end), sum(ROLLING_WEIGHTS_CONST_RETURN(:, 1:7), 2) * 100, 'DisplayName', 'Equity');
hold on;
plot(DATES(WINDOW_SIZE:end), sum(ROLLING_WEIGHTS_CONST_RETURN(:, 8:10), 2) * 100, 'DisplayName', 'Bond');
hold off;
title('Portfolio composition (constant return, variable risk, annualized)');
xlabel('Year');
ylabel('Weight (%)');
ylim([0 100]); % Limita l'asse Y a 100
legend('show');
%% Focus sugli ultimi 5 anni (dal 2018 al 2023) - Seconda Analisi Modificata
if ~isempty(focus_indices)
    FOCUS_RISKS_CONST_RETURN = ROLLING_RISKS_CONST_RETURN(focus_indices);
    FOCUS_WEIGHTS_CONST_RETURN = ROLLING_WEIGHTS_CONST_RETURN(focus_indices, :);

    figure;
    subplot(2,1,1);
    plot(FOCUS_DATES, FOCUS_RISKS_CONST_RETURN * 100);
    title('Risk (2018-2023, variable)');
    xlabel('Year');
    ylabel('Risk (%)');

    subplot(2,1,2);
    plot(FOCUS_DATES, sum(FOCUS_WEIGHTS_CONST_RETURN(:, 1:7), 2) * 100, 'DisplayName', 'Equity');
    hold on;
    plot(FOCUS_DATES, sum(FOCUS_WEIGHTS_CONST_RETURN(:, 8:10), 2) * 100, 'DisplayName', 'Bond');
    hold off;
    title('Portfolio composition (2018-2023, variable risk)');
    xlabel('Year');
    ylabel('Weight (%)');
    ylim([0,100])
    legend('show');
end
