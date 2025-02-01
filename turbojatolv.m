%% Prova 3 - Propulsão 1
% Cálculo de motor turbojato ideal ou não-ideal

%% Inputs
ideal       = 0;                % 1 para caso ideal e 0 para não ideal
M_0         = 0.8;                % Mach de voo
h           = 35000*0.3048;    % Altitude de voo (caso fornecida) --> Converteru de pés para metros
    % [T_0, ~, P_0, ~] = atmosisa(h);
    T_0 = 288 ;
    P_0 = 101.3 ; 
    % P_0 = P_0/1000; %[kPa]

gamma_c    = 1.4;      % Razao dos calores especificos do ar no compressor
cp_c        = 1.005;    % Calor especifico no compressor [kJ/kg.K]
h_PR        = 44000;    % Poder calorífico do combustível [kJ/kg]
pi_dmax     = 0.99;     % Razão de pressão no difusor


if ideal
    pi_b        = 1;
    pi_n        = 1;
    eta_n = pi_n;
    e_c         = 1;
    e_t         = 1;
    eta_b       = 1;        
    eta_m       = 1;
    gamma_t    = gamma_c;
    cp_t        = cp_c;
   
else
    pi_b        = 0.98;         % Razão de pressão na camera de combustão
    pi_n        = 1;         % Razão de pressão no bocal de saida
    eta_n = pi_n;
    eta_d = 0.9;
    eta_c = 0.89; % Eficiencia total de compressão
    e_c         = 1;          % Eficiencia politrópica do compressor
    e_t         = 1;          % Eficiencia politrópica da turbina
    eta_b       = 0.98;         % Eficiencia da combustão
    eta_m       = 0.95;         % Eficiencia mecanica
    gamma_t    = 1.3;          % Razao dos calores especificos do ar na turbina
    cp_t        = 1.239;        % Calor especifico na turbina
    eta_t = 0.98;
end

% Variáveis de controle

pi_c        = 25;       % Razão de pressão no compressor    (aumentar diminui o consumo)
Tt_4        = 2000;     % Temperatura na entrada da turbina (diminuir diminui o consumo)

g_c         = 1;        % Sistema imperial

% % Análise de chocked/unchocked
% P0_P9       = 0.9;      % Razão de pressão no bocal (P0/P9)
% P0_P19      = 0.9;      % Razão de pressão no bocal do fan (P0/P19)


%% Outputs
% F/dot_m0  - Empuxo específico
% f         - Razão entre massa de combustível e massa de ar
% S         - Consumo específico
% eta_T     - Eficiencia térmica
% eta_P     - Eficiencia propulsiva
% eta_0     - Eficiencia total
% eta_c     - Eficiencia do compressor
% eta_t     - Eficiencia da turbina
%

%% Equations
% Definições iniciais
% 7-52a
R_c     = ((gamma_c-1)/gamma_c) * cp_c;
% 7-52b
R_t     = ((gamma_t-1)/gamma_t) * cp_t;
% 7-52c
a_0     = sqrt(gamma_c*(R_c*1000)*g_c*T_0);
% 7-52d
V_0     = a_0*M_0;

% Far away (r)
% 7-52e
tal_r   = 1 + (gamma_c-1)/2 * M_0^2;
% 7-52f
pi_r    = tal_r^(gamma_c/(gamma_c-1));

if M_0 > 1
    fprintf('Mach de entrada é maior que 1\n')
    % 7-52h
    eta_r   = 1 - 0.075*(M_0-1)^1.35;
else
    % 7-52g
    eta_r   = 1;
end

% Difusor (d)
% 7-52i (razão de pressão no difusor (Pt2/Pt1)
if ideal
    pi_d = 1;
else
    pi_d = pi_dmax*eta_d;
end

% tal_d = pi_d^((gamma_c-1)/gamma_c);


% Compressor (c)
pi_c = pi_c * eta_c;

% 7-52k
tal_c       = pi_c^((gamma_c-1)/(gamma_c*e_c));
% 7-52l
% eta_c       = (pi_c^((gamma_c-1)/gamma_c) - 1) / (tal_c-1);

% Combustor (b)
% 7-52j
tal_lambda  = cp_t*Tt_4/(cp_c*T_0);
% 7-52o
f           = (tal_lambda - tal_r*tal_c) / ((h_PR*eta_b / (cp_c*T_0)) - tal_lambda);

% Turbine (t)
% 7-52p
tal_t       = eta_t*(1 - ((1/(eta_m*(1+f))) * (tal_r/tal_lambda) * (tal_c - 1 )));
% 7-52q
pi_t        = tal_t^(gamma_t/((gamma_t-1)*e_t));
% 7-52r
% eta_t       = (1-tal_t) / (1-tal_t^(1/e_t));

%% Chocked or Unchocked
if 1
    Pt_5 = pi_r*pi_d*pi_c*pi_b*pi_t*pi_n*P_0;
    Pc_C = Pt_5*((1-(1/eta_n)*((gamma_t-1)/(gamma_t+1)))^(gamma_t/(gamma_t-1)));

    % P0_P19 = (((gamma_c+1)/2)^(gamma_c/(gamma_c-1))) / (pi_r*pi_d*pi_f*pi_fn);
    % P0_P9  = (((gamma_t+1)/2)^(gamma_t/(gamma_t-1))) / (pi_r*pi_d*pi_c*pi_b*pi_t*pi_n);

    if Pc_C > P_0 
        fprintf('Nozzle Choked!\n')
        P0_P9   = P_0/Pc_C;
        

        % Core
        Tt_5 = tal_t*Tt_4;
        T_9 = Tt_5/((gamma_t+1)/2);
        V_9 = sqrt(gamma_t*R_t*1000*T_9);
        T9_T0 = T_9/T_0;
        V9_a0 = V_9/a_0;
        

    else
        fprintf('Nozzle Unchoked!\n')
        P0_P9  = 1;
       
        % 7-52s
        Pt9_P9      = (P0_P9) * pi_r*pi_d*pi_c*pi_b*pi_t*pi_n;
        % 7-52t
        M_9         = sqrt((2/(gamma_t-1)) * ((Pt9_P9)^((gamma_t-1)/gamma_t) - 1));
        % 7-52u
        T9_T0       = ((tal_lambda*tal_t)/((Pt9_P9)^((gamma_t-1)/gamma_t)))*(cp_c/cp_t);
        % 7-52v
        V9_a0       = M_9 * sqrt(((gamma_t*R_t)/(gamma_t*R_c))*T9_T0);
        
    end
end
%%
F_dot_m0_core = (a_0/g_c) * ((1+f)*V9_a0 - M_0 + (1+f) *(((R_t*T9_T0)/(R_c*V9_a0))*((1-(P0_P9))/gamma_c)));


% 7-52aa
F_dot_m0    = (a_0/g_c) * ((1+f)* V9_a0 - M_0 + (1+f) * (((R_t*T9_T0)/(R_c* V9_a0)) * ((1-(P0_P9))/gamma_c))) ;
% 7-52ab
S           = f/F_dot_m0; %[kg/s / N]

% 7-52ad
eta_T       = (a_0^2 * ((1 + f) * (V9_a0)^2 -  M_0^2)) / (2 * g_c * f * h_PR * 1000); % Eficiência térmica
% 7-52ae
eta_P = (2 * g_c * V_0 * F_dot_m0)/(a_0^2 * ( (1+f)* (V9_a0^2) - M_0^2)); % Eficiência propulsiva
% 7-52af
eta_0       = eta_P*eta_T;


% Checks adicionais
V_9 = V9_a0*a_0;
T_9 = T9_T0*T_0;

Tt_5 = tal_t*Tt_4;
T_9y = Tt_5/((gamma_t+1)/2);
V_9y = sqrt(gamma_t*R_t*1000*T_9y);
F = F_dot_m0*755/1000; % [kN]

% minimum specific fuel consumption and maximum specific thrust
% tal_f_star = (tal_lambda * tal_r * (tal_c - 1) - tal_lambda / (tal_r * tal_c) + alfa * tal_r + 1) / ...
%              (tal_r * (1 + alfa));
% 
% pi_f_star = tal_f_star^(gamma_c*e_f/(gamma_c-1));