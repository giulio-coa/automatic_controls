%% Esercizio #1
close all
clear all

% Dati
s = tf('s');

F_1 = (1 + 10 * s) / ((1 + 0.1 * s) * (1 + 5 * s));
F_2 = 1 / s;
K_r = 1;
d = 1.5;

% Analisi delle specifiche statiche
K_F1 = dcgain(F_1);
K_F2 = dcgain(s * F_2);

h = 1
K_c = 6.25

% Analisi delle specifiche dinamiche
w_bdes = 4;
w_cdes = 0.63 * w_bdes
w_c = 2.2;

s_max = 0.25;
M_r = (1 + s_max) / 0.85;
M_rdB = 20 * log10(M_r);
m_phi = 60 - 5 * M_rdB
m_phi = 50;

% Prima approssimazione della funzione d'anello
G_a1 = (K_c * F_1 * F_2) / (s^h * K_r);
[module_1, phase_1] = bode(G_a1, w_c);
recovery = m_phi - phase_1 - 180;

% Rete derivatrice
m_d1 = 3;
x_d1 = sqrt(m_d1);
tau_d1 = x_d1 / w_c;
R_d1 = (1 + tau_d1 * s) / (1 + tau_d1 / m_d1 * s);

m_d2 = 4;
x_d2 = sqrt(m_d2);
tau_d2 = x_d2 / w_c;
R_d2 = (1 + tau_d2 * s) / (1 + tau_d2 / m_d2 * s);

R_d = R_d1 * R_d2;

G_a2 = G_a1 * R_d;
[module_2, phase_2] = bode(G_a2, w_c);

% Rete integratrice
m_i = 9;
x_i = 150;
tau_i = x_i / w_c;
R_i = (1 + tau_i / m_i * s) / (1 + tau_i * s);

G_a3 = G_a2 * R_i;
[module_3, phase_3] = bode(G_a3, w_c);

% Controllore
C = (K_c * R_d * R_i) / s^h;

G_a = C * F_1 * F_2 / K_r;
K_Ga = dcgain(s^(h + 1) * G_a);

figure, margin(G_a)

% Verifica delle specifiche
e_r = abs(K_r / K_Ga)

W = feedback(C * F_1 * F_2, 1 / K_r);

figure, bode(W / dcgain(W))
figure, step(W / dcgain(W))

% Valutazioni extra
u = 1.052 % runna schema Simulink e trova il valore massimo dallo scope di u

% Discretizzazione
w_b = 4.28;
T_1 = (2 * pi) / (20 * w_b)
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

F = (100 * (s + 10)) /(s * (s^4 + 38 * s^3 + 481 * s^2 + 2280 * s + 3600));
N = 10;

[m_G, m_phi, w_mG, w_mphi] = margin(F);
w_pi = 3.5357;
[module, phase] = bode(F, w_pi)

% Calcolo dei parametri per il metodo Ziegler-Nichols in catena chiusa
K_p_ = m_G;
T_ = 2 * pi / w_pi;

% Calcolo dei parametri per il regolatore PID
K_P = 0.6 * K_p_
T_i = 0.5 * T_
T_d = 0.125 * T_

R = K_P * (1 + 1 / (T_i * s) + (T_d * s) / (1 + T_d / N * s));

W = feedback(R * F, 1);
figure, bode(W / dcgain(W))
