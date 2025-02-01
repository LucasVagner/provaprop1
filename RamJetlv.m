%% Prova 3 - Propulsão I - CÁLCULO DE MOTOR RAMJET
clear; clc;

%% Dados de Entrada
M_0 = 2.4;
T_0 = 216.7; % [K]
gamma = 1.4;
cp = 1004; % [J/kg.K]
h_PR = 42800000; % [J/kg]
Tt_4 = 2200; % [K]
g_c = 1;

%% Cálculos
% Calcular R
R = ((gamma-1)/gamma) * cp;

% Calcular a0
a_0 = sqrt(gamma*R*g_c*T_0);

% Calcular tau_r
tau_r = 1 + ((gamma-1)/2)*M_0^2;

% Calcular tau_lambda
tau_lambda = Tt_4/T_0;

% Calcular V9/a0
V9_a0 = M_0 * sqrt(tau_lambda/tau_r);

% Calcular F/m0 - Empuxo Específico
F_m0 = (a_0/g_c) * (V9_a0 - M_0);

% Calcular f
f = (cp* T_0/h_PR) * (tau_lambda-tau_r);

% Calcular S - Consumo Específico
S = f/F_m0;

% Calcular eta_T - Eficiência Térmica
eta_T = 1 - (1/tau_r);

% Calcular eta_P - Eficiência Propulsiva
eta_P = 2 / (sqrt(tau_lambda/tau_r) + 1);

% Calcular eta_O - Eficiência Total
eta_O = eta_T*eta_P;

%% Display Results in Table Format
resultados = table(...
    ["R"; "a_0"; "tau_r"; "tau_lambda"; "V9_a0"; "F_m0"; "f"; "S"; "eta_T"; "eta_P"; "eta_O"], ...
    ["Constante do ar"; "Velocidade do som na entrada"; "Razão de temperatura de estagnação a montante"; ...
     "Razão de temperatura de estagnação na turbina"; "Razão de velocidade de saída"; "Empuxo específico"; ...
     "Fração de combustível"; "Consumo específico"; "Eficiência térmica"; "Eficiência propulsiva"; "Eficiência total"], ...
    [R; a_0; tau_r; tau_lambda; V9_a0; F_m0; f; S; eta_T; eta_P; eta_O], ...
    ["J/(kg*K)"; "m/s"; "-"; "-"; "-"; "N/(kg/s)"; "-"; "kg/(N*s)"; "-"; "-"; "-"], ...
    'VariableNames', {'Símbolo', 'Variável', 'Valor', 'Unidade'});

% Exibição do DataFrame final
disp(resultados);