%% Prova 3 - Propulsão 1
% Cálculo de motor turbofan ideal ou não-ideal

%% Inputs
ideal       = 0;                % 1 para caso ideal e 0 para não ideal
M_0         = 0.8;                % Mach de voo
h           = 35000*0.3048;    % Altitude de voo (caso fornecida) --> Converteru de pés para metros
    % [T_0, ~, P_0, ~] = atmosisa(h);
    T_0 = 288 ;
    P_0 = 101.3 ; 
    % P_0 = P_0/1000; %[kPa]

gamma_c    = 1.4;      % Razao dos calores especificos do ar no compressor
cp_c        = 1.004;    % Calor especifico no compressor [kJ/kg.K]
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
    pi_fn       = 1;      % fan nozzle total pressure ratio Pt19/Pt13
    eta_fn = 1;
    e_f         = 1; 

else
    pi_b        = 0.98;         % Razão de pressão na camera de combustão
    pi_n        = 1;         % Razão de pressão no bocal de saida
    eta_n = pi_n;
    pi_fn       = 1;         % Fan nozzle total pressure ratio Pt19/Pt13
    eta_fn = pi_fn;
    eta_d = 0.9;
    eta_c = 0.89;            % Eficiencia total de compressão
    e_f         = 1;         % Eficiencia politropica do fan
    e_c         = 1;         % Eficiencia politrópica do compressor
    e_t         = 1;         % Eficiencia politrópica da turbina
    eta_b       = 0.98;      % Eficiencia da combustão
    eta_m       = 0.95;        % Eficiencia mecanica
    eta_f       = 0.91;         % Eficiencia do fan
    gamma_t    = 1.3;          % Razao dos calores especificos do ar na turbina
    cp_t        = 1.239;        % Calor especifico na turbina
    eta_t = 0.98;
end

% Variáveis de controle

pi_f        = 3.375;      % Fan total pressure ratio Pt13/Pt2 (aumentar diminui o consumo)
pi_c        = 25;       % Razão de pressão no compressor    (aumentar diminui o consumo)
pi_core     = pi_c/pi_f;  % Razão de pressão no core    
alfa        = 0.7;       % Razão de bypass                   (aumentar diminui o consumo)
Tt_4        = 2000;     % Temperatura na entrada da turbina (diminuir diminui o consumo)

g_c         = 1;        % Sistema imperial

% Análise de chocked/unchocked
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
tal_d = pi_d^((gamma_c-1)/gamma_c);

% Fan (f)
pi_f = pi_f*eta_f;
% 7-52m
tal_f       = pi_f^((gamma_c-1)/(gamma_c*e_f));
% 7-52n
% eta_f       = (pi_f^((gamma_c-1)/gamma_c) - 1) / (tal_f-1);


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
tal_t       = eta_t*(1 - ((1/(eta_m*(1+f)))*(tal_r/tal_lambda)*(tal_c - 1 + alfa*(tal_f-1))));
% 7-52q
pi_t        = tal_t^(gamma_t/((gamma_t-1)*e_t));
% 7-52r
% eta_t       = (1-tal_t) / (1-tal_t^(1/e_t));

%% Chocked or Unchocked
if 1
    Pt_5 = pi_r*pi_d*pi_c*pi_b*pi_t*P_0;
    Pc_C = Pt_5*((1-(1/eta_n)*((gamma_t-1)/(gamma_t+1)))^(gamma_t/(gamma_t-1)));

    Pt_13 = pi_r*pi_d*pi_f*P_0;
    Pc_F  = Pt_13*((1-(1/eta_fn)*((gamma_c-1)/(gamma_c+1)))^(gamma_c/(gamma_c-1))); 

    % P0_P19 = (((gamma_c+1)/2)^(gamma_c/(gamma_c-1))) / (pi_r*pi_d*pi_f*pi_fn);
    % P0_P9  = (((gamma_t+1)/2)^(gamma_t/(gamma_t-1))) / (pi_r*pi_d*pi_c*pi_b*pi_t*pi_n);

    if Pc_C > P_0 || Pc_F > P_0
        fprintf('Nozzle Choked!\n')
        P0_P9   = P_0/Pc_C;
        P0_P19  = P_0/Pc_F;

        % Core
        Tt_5 = tal_t*Tt_4;
        T_9 = Tt_5/((gamma_t+1)/2);
        V_9 = sqrt(gamma_t*R_t*1000*T_9);
        T9_T0 = T_9/T_0;
        V9_a0 = V_9/a_0;
        
        % Fan
        Tt_13 = tal_d*tal_r*tal_f*T_0;
        T_19 = Tt_13/((gamma_c+1)/2);
        V_19 = sqrt(gamma_c*R_c*1000*T_19);
        T19_T0 = T_19/T_0;
        V19_a0 = V_19/a_0;

    else
        fprintf('Nozzle Unchoked!\n')
        P0_P9  = 1;
        P0_P19 = 1;
        % 7-52s
        Pt9_P9      = (P0_P9) * pi_r*pi_d*pi_c*pi_b*pi_t*pi_n;
        % 7-52t
        M_9         = sqrt((2/(gamma_t-1)) * ((Pt9_P9)^((gamma_t-1)/gamma_t) - 1));
        % 7-52u
        T9_T0       = ((tal_lambda*tal_t)/((Pt9_P9)^((gamma_t-1)/gamma_t)))*(cp_c/cp_t);
        % 7-52v
        V9_a0       = M_9 * sqrt(((gamma_t*R_t)/(gamma_t*R_c))*T9_T0);
        % 7-52w
        Pt19_P19    = (P0_P19) * pi_r*pi_d*pi_f*pi_fn;
        % 7-52x
        M_19        = sqrt((2/(gamma_c-1)) * ((Pt19_P19)^((gamma_c-1)/gamma_c) - 1));
        % 7-52y
        T19_T0      = (tal_r*tal_f)/((Pt19_P19)^((gamma_c-1)/gamma_c));
        % 7-52z
        V19_a0      = M_19 * sqrt(T19_T0);
    end
end
%%
F_dot_m0_core = (1/(alfa+1)) * (a_0/g_c) * ((1+f) * V9_a0 - M_0 + (1+f) *(((R_t*T9_T0)/(R_c*V9_a0))*((1-(P0_P9))/gamma_c)));

F_dot_m0_fan = (alfa/(alfa + 1)) * (a_0/g_c) * (V19_a0 - M_0 + (T19_T0/V19_a0)*((1-(P0_P19))/gamma_c));

F_dot_m0_X = F_dot_m0_fan + F_dot_m0_core ;

% 7-52aa
F_dot_m0    = (1/(1+alfa)) * (a_0/g_c) * ((1+f)*V9_a0 - M_0 + (1+f) *(((R_t*T9_T0)/(R_c*V9_a0))*((1-(P0_P9))/gamma_c))) + (alfa/(1+alfa))*(a_0/g_c)...
    *(V19_a0 - M_0 + (T19_T0/V19_a0)*((1-(P0_P19))/gamma_c)); %[N / kg/s]

% 7-52ab
S           = f/((1+alfa)*(F_dot_m0)); %[kg/s / N]
% 7-52ac
FR = ((1+f)*(V9_a0 - M_0) + (1+f) * (R_t/R_c) * (T9_T0/V9_a0)*((1-P0_P9)/gamma_c)) / ((V19_a0 - M_0) + (T19_T0/V19_a0)*((1-P0_P19)/gamma_c));

% 7-52ad
eta_T       = (a_0^2 * ((1 + f) * (V9_a0)^2 + alfa * (V19_a0)^2 - (1 + alfa) * M_0^2)) / (2 * g_c * f * h_PR * 1000);
% 7-52ae
eta_P = (2 * M_0 * ((1 + f) * (V9_a0) + alfa * (V19_a0) - (1 + alfa) * M_0)) / ...
        ((1 + f) * (V9_a0)^2 + alfa * (V19_a0)^2 - (1 + alfa) * M_0^2);
% 7-52af
eta_0       = eta_P*eta_T;


% Checks adicionais
V_9 = V9_a0*a_0;
T_9 = T9_T0*T_0;

Tt_5 = tal_t*Tt_4;
T_9y = Tt_5/((gamma_t+1)/2);
V_9y = sqrt(gamma_t*R_t*1000*T_9y);
F = F_dot_m0*755/1000; % [kN]

% bypass ratio for the minimum fuel consumption
alfa_star = (1 / (tal_r * (tal_f - 1))) * ...
    ((tal_lambda - tal_r * (tal_c - 1)) - (tal_lambda / (tal_r * tal_c)) - ...
    (1/4) * (sqrt(tal_r * tal_f - 1) + sqrt(tal_r - 1))^2);

% minimum specific fuel consumption and maximum specific thrust
tal_f_star = (tal_lambda * tal_r * (tal_c - 1) - tal_lambda / (tal_r * tal_c) + alfa * tal_r + 1) / ...
             (tal_r * (1 + alfa));

pi_f_star = tal_f_star^(gamma_c*e_f/(gamma_c-1));