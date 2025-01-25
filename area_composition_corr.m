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
WINDOW_SIZE = 3 * 12; % Finestra temporale di 3 anni (in mesi)
NUM_WINDOWS = size(RETURNS, 1) - WINDOW_SIZE + 1; % Numero di finestre temporali

% Inizializza matrici per registrare i risultati
ROLLING_WEIGHTS = zeros(NUM_WINDOWS, ASSET);
ROLLING_RETURNS = zeros(NUM_WINDOWS, 1);
ROLLING_RISKS = zeros(NUM_WINDOWS, 1);

% Simulazione Monte Carlo con resampling per ogni finestra temporale
for w = 1:NUM_WINDOWS
    % Estrai i dati per la finestra temporale corrente
    WINDOW_RETURNS = RETURNS(w:w + WINDOW_SIZE - 1, :);
    
    % Controllo NaN
    if any(isnan(WINDOW_RETURNS), 'all')
        disp(['Valori NaN trovati nella finestra ', num2str(w)]);
        continue; % Salta a prossima iterazione
    end

    % Calcola i ritorni attesi annuali e la matrice di covarianza annuale per la finestra temporale corrente
    mean_returns = mean(WINDOW_RETURNS, 'omitnan'); % 'omitnan' per ignorare NaN
    EXP_RET_annual = (1 + mean_returns).^12 - 1; % Ritorni annualizzati
    COV_annual = cov(WINDOW_RETURNS, 'omitrows') * 12; % Matrice di covarianza annualizzata

    % Ottieni i pesi del portafoglio ottimizzato con i dati della finestra corrente
    [RISK, ROR, WTS] = portopt(EXP_RET_annual', COV_annual, 1); % Ottimizzazione per 1 portafoglio

    % Memorizza i risultati della finestra corrente
    if size(WTS, 2) == ASSET
        ROLLING_WEIGHTS(w, :) = WTS; % Assicurati che le dimensioni siano corrette
    else
        error('Dimensione di WTS non corrisponde al numero di asset');
    end
    
    ROLLING_RETURNS(w) = ROR;
    ROLLING_RISKS(w) = RISK;
end

% Verifica le dimensioni delle date e dei dati da plottare
disp('Dimensioni delle date:');
disp(size(DATES(WINDOW_SIZE:end)));
disp('Dimensioni dei ritorni rolling:');
disp(size(ROLLING_RETURNS));
disp('Dimensioni dei rischi rolling:');
disp(size(ROLLING_RISKS));
disp('Dimensioni dei pesi rolling:');
disp(size(ROLLING_WEIGHTS));

% Calcola i pesi dei portafogli ottimali tramite resampling
RESAPL_WEIGHTS = mean(ROLLING_WEIGHTS, 1);

% Imposta la percentuale desiderata
percent_equity = 0.6; % 60% per gli asset di equity
percent_bond = 0.4; % 40% per i bond

% Normalizza i pesi per rispettare la percentuale desiderata
total_equity_weight = sum(RESAPL_WEIGHTS(1:7));
RESAPL_WEIGHTS(1:7) = RESAPL_WEIGHTS(1:7) * (percent_equity / total_equity_weight);

total_bond_weight = sum(RESAPL_WEIGHTS(8:10));
RESAPL_WEIGHTS(8:10) = RESAPL_WEIGHTS(8:10) * (percent_bond / total_bond_weight);

% Verifica che la somma sia corretta
check_equity = sum(RESAPL_WEIGHTS(1:7)); % Somma dei pesi degli asset di equity per ogni simulazione
check_bond = sum(RESAPL_WEIGHTS(8:10)); % Somma dei pesi dei bond per ogni simulazione

disp(['Media della somma dei pesi di equity: ', num2str(check_equity)]);
disp(['Media della somma dei pesi di bond: ', num2str(check_bond)]);

avg_weights = RESAPL_WEIGHTS; % Media lungo le righe, quindi ottieni la media per ogni colonna (asset)

disp('Average weights of the portfolio:');
disp(['Equity: ', num2str(avg_weights(1:7))]);
disp(['Bond: ', num2str(avg_weights(8:10))]);

% Grafici
figure;
subplot(3,1,1);
plot(DATES(WINDOW_SIZE:end), ROLLING_RETURNS * 100);
title('Expected return of the portfolio in time');
xlabel('Year');
ylabel('Expected Return (%)');

subplot(3,1,2);
plot(DATES(WINDOW_SIZE:end), ROLLING_RISKS * 100);
title('Risk of the portfolio in time');
xlabel('Year');
ylabel('Risk (%)');

subplot(3,1,3);
plot(DATES(WINDOW_SIZE:end), sum(ROLLING_WEIGHTS(:, 1:7), 2) * 100, 'DisplayName', 'Equity');
hold on;
plot(DATES(WINDOW_SIZE:end), sum(ROLLING_WEIGHTS(:, 8:10), 2) * 100, 'DisplayName', 'Bond');
hold off;
title('Portfolio composition in time');
xlabel('Year');
ylabel('Weight (%)');
legend('show');
%%
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
WINDOW_SIZE = 3 * 12; % Finestra temporale di 3 anni (in mesi)
NUM_WINDOWS = size(RETURNS, 1) - WINDOW_SIZE + 1; % Numero di finestre temporali

% Inizializza matrici per registrare i risultati
ROLLING_WEIGHTS = zeros(NUM_WINDOWS, ASSET);
ROLLING_RETURNS = zeros(NUM_WINDOWS, 1);
ROLLING_RISKS = zeros(NUM_WINDOWS, 1);

% Simulazione Monte Carlo con resampling per ogni finestra temporale (cambiando rendimenti e correlazioni)
for w = 1:NUM_WINDOWS
    % Estrai i dati per la finestra temporale corrente
    WINDOW_RETURNS = RETURNS(w:w + WINDOW_SIZE - 1, :);
    
    % Calcola i ritorni attesi annuali e la matrice di covarianza annuale per la finestra temporale corrente
    mean_returns = mean(WINDOW_RETURNS);
    EXP_RET_annual = (1 + mean_returns).^12 - 1; % Ritorni annualizzati
    COV_annual = cov(WINDOW_RETURNS) * 12; % Matrice di covarianza annualizzata
    
    % Ottieni i pesi del portafoglio ottimizzato
    [RISK, ROR, WTS] = portopt(EXP_RET_annual', COV_annual, 1); % Ottimizzazione per 1 portafoglio
    
    % Memorizza i risultati della finestra corrente
    ROLLING_WEIGHTS(w, :) = WTS;
    ROLLING_RETURNS(w) = ROR;
    ROLLING_RISKS(w) = RISK;
end

% Calcola i pesi medi e rischi per la seconda analisi (tenendo fissi i rendimenti e i rischi)
% Calcola i valori medi da utilizzare
mean_returns_fixed = mean(ROLLING_RETURNS); % Rendimento medio
mean_risk_fixed = mean(ROLLING_RISKS); % Rischio medio

% Inizializza matrici per i risultati della seconda analisi
ROLLING_WEIGHTS_FIXED = zeros(NUM_WINDOWS, ASSET);
ROLLING_RETURNS_FIXED = zeros(NUM_WINDOWS, 1);
ROLLING_RISKS_FIXED = zeros(NUM_WINDOWS, 1);

% Simulazione Monte Carlo con resampling per ogni finestra temporale (solo correlazioni)
for w = 1:NUM_WINDOWS
    % Estrai i dati per la finestra temporale corrente
    WINDOW_RETURNS = RETURNS(w:w + WINDOW_SIZE - 1, :);
    
    % Calcola la matrice di covarianza annuale per la finestra temporale corrente
    COV_annual_fixed = cov(WINDOW_RETURNS) * 12; % Matrice di covarianza annualizzata

    % Mantieni fissi i rendimenti e rischi medi
    [RISK_FIXED, ROR_FIXED, WTS_FIXED] = portopt(mean_returns_fixed', COV_annual_fixed, 1); % Ottimizzazione per 1 portafoglio
    
    % Memorizza i risultati della finestra corrente
    ROLLING_WEIGHTS_FIXED(w, :) = WTS_FIXED;
    ROLLING_RETURNS_FIXED(w) = ROR_FIXED;
    ROLLING_RISKS_FIXED(w) = RISK_FIXED;
end

% Grafici per entrambe le analisi
figure;

% Grafico 1: Ritorno nel tempo - Tutto variabile
subplot(2,2,1);
plot(DATES(WINDOW_SIZE:end), ROLLING_RETURNS * 100, 'r');
title('Rendimento del portafoglio - Tutto variabile');
xlabel('Anno');
ylabel('Rendimento (%)');

% Grafico 2: Rischio nel tempo - Tutto variabile
subplot(2,2,2);
plot(DATES(WINDOW_SIZE:end), ROLLING_RISKS * 100, 'r');
title('Rischio del portafoglio - Tutto variabile');
xlabel('Anno');
ylabel('Rischio (%)');

% Grafico 3: Ritorno nel tempo - Solo correlazioni
subplot(2,2,3);
plot(DATES(WINDOW_SIZE:end), ROLLING_RETURNS_FIXED * 100, 'b');
title('Rendimento del portafoglio - Solo correlazioni');
xlabel('Anno');
ylabel('Rendimento (%)');

% Grafico 4: Rischio nel tempo - Solo correlazioni
subplot(2,2,4);
plot(DATES(WINDOW_SIZE:end), ROLLING_RISKS_FIXED * 100, 'b');
title('Rischio del portafoglio - Solo correlazioni');
xlabel('Anno');
ylabel('Rischio (%)');

% Grafico della composizione dei pesi per ciascuna simulazione
figure;

% Grafico composizione - Tutto variabile
subplot(2,1,1);
area(DATES(WINDOW_SIZE:end), ROLLING_WEIGHTS);
title('Composizione del portafoglio - Tutto variabile');
xlabel('Anno');
ylabel('Pesi (%)');
legend(LABELS, 'Location', 'EastOutside');

% Grafico composizione - Solo correlazioni
subplot(2,1,2);
area(DATES(WINDOW_SIZE:end), ROLLING_WEIGHTS_FIXED);
title('Composizione del portafoglio - Solo correlazioni');
xlabel('Anno');
ylabel('Pesi (%)');
legend(LABELS, 'Location', 'EastOutside');
