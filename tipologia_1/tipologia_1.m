%% Esercizio #1
close all
clear all

% Dati
s = tf('s');

F_1 = 30 / (s + 15);
F_2 = (3 * (s + 1)) / (s * (s + 4) * (s + 6));
K_r = 1;
d_1 = 1;
d_2 = 4;

% Analisi delle specifiche statiche
K_F1 = dcgain(F_1);
K_F2 = dcgain(s * F_2);

h = 0
K_c = 40

% Analisi delle specifiche dinamiche
w_b = 20;
w_cdes = 0.63 * w_b
w_c = 13;

s_max = 0.2;
M_r = (1 + s_max) / 0.85;
M_rdB = 20 * log10(M_r);
m_phi = 60 - 5 * M_rdB

% Prima approssimazione della funzione d'anello
G_a1 = (K_c * F_1 * F_2) / (s^h * K_r);
[module_1, phase_1] = bode(G_a1, w_c);
recovery = m_phi - phase_1 - 180;

% Rete derivatrice
m_d = 12;
x_d = sqrt(m_d);
tau_d = x_d / w_c;
R_d = (1 + tau_d * s) / (1 + tau_d / m_d * s);

G_a2 = G_a1 * R_d;
[module_2, phase_2] = bode(G_a2, w_c);

% Rete integratrice
m_i = 3.5;
x_i = 50;
tau_i = x_i / w_c;
R_i = (1 + tau_i / m_i * s) / (1 + tau_i * s);

G_a3 = G_a2 * R_i;
[module_3, phase_3] = bode(G_a3, w_c);

% Controllore
C = (K_c * R_d * R_i) / s^h

G_a = C * F_1 * F_2 / K_r;
K_Ga = dcgain(s^(h + 1) * G_a);

figure, margin(G_a)

% Verifica delle specifiche
e_r = abs(K_r / K_Ga)

e_d1 = abs(d_1 / dcgain(C * F_1))

e_d2 = abs(d_2 / K_Ga)

W = feedback(C * F_1 * F_2, 1 / K_r);

figure, step(W / dcgain(W))
figure, bode(W / dcgain(W))

% Valutazioni extra
w_0 = 0.2;
amplitude = 1;
W_e = feedback(1, C * F_1 * F_2 / K_r);
[module_e, phase_e] = bode(W_e, w_0);
e_max = module_e * amplitude

% Discretizzazione
w_b = 21;
T_1 = (2 * pi) / (20 * w_b);
T = T_1
G_aZOH = G_a / (1 + T / 2 * s);
C_z = c2d(C, T, 'tustin')

figure, margin(G_aZOH)

F_1z = c2d(F_1, T, 'tustin');
F_2z = c2d(F_2, T, 'tustin');
W_z = feedback(C_z * F_1z * F_2z, 1 / K_r);

figure, step(W_z / dcgain(W_z))

%% Esercizio #3
close all
clear all

% Dati
s = tf('s');

F = (4 * (s + 150)^2) /((s + 100) * (s + 4) * (s + 50));
N = 10;
y_inf = 4.5;
u = 1;

% Calcolo dei parametri per il metodo Ziegler-Nichols in catena aperta
K_F = y_inf / u
theta_F = 0.01
y_theta_F_tau_F = 0.63 * y_inf;
tau_F = 0.26 - theta_F

% Calcolo dei parametri per il regolatore PID
K_P = (1.2 * tau_F) / (K_F * theta_F)
T_i = 2 * theta_F
T_d = 0.5 * theta_F

R = K_P * (1 + 1 / (T_i * s) + (T_d * s) / (1 + T_d / N * s));

W = feedback(R * F, 1);
figure, step(W / dcgain(W))
