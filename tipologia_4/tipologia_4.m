%% Esercizio #1
close all
clear all

% Dati
s = tf('s');

F_1 = 5 / s;
F_2 = (s + 20) / ((s + 1) * (s + 5)^2);
K_r = 1;
d_1 = 0.5;
d_2 = 0.1;

% Analisi delle specifiche statiche
K_F1 = dcgain(s * F_1);
K_F2 = dcgain(F_2);

h = 0
K_c = 5

% Analisi delle specifiche dinamiche
t_s = 1;
w_bdes = 3 / t_s;
w_cdes = 0.63 * w_bdes
w_c = 2;

M_rdB = 2.5;
m_phi = 60 - 5 * M_rdB
m_phi = 50;

% Prima approssimazione della funzione d'anello
G_a1 = (K_c * F_1 * F_2) / (s^h * K_r);
[module_1, phase_1] = bode(G_a1, w_c);
recovery = m_phi - phase_1 - 180;

% Rete derivatrice
m_d = 4;
x_d = sqrt(m_d);
tau_d = x_d / w_c;
R_d = (1 + tau_d * s) / (1 + tau_d / m_d * s);

R_d = R_d^2;

G_a2 = G_a1 * R_d;
[module_2, phase_2] = bode(G_a2, w_c);

% Rete integratrice
m_i = 10;
x_i = 100;
tau_i = x_i / w_c;
R_i = (1 + tau_i / m_i * s) / (1 + tau_i * s);

R_i = R_i^2;

G_a3 = G_a2 * R_i;
[module_3, phase_3] = bode(G_a3, w_c);

% Fattore correttivo
K_correttivo = 1 / module_3

K_c = K_c * K_correttivo;

G_a4 = G_a3 * K_correttivo;
[module_4, phase_4] = bode(G_a4, w_c);

% Controllore
C = (K_c * R_d * R_i) / s^h;

G_a = C * F_1 * F_2 / K_r;
K_Ga = dcgain(s^(h + 1) * G_a);

figure, margin(G_a)

% Verifica delle specifiche
e_r = abs(K_r / K_Ga)

e_d2 = abs(d_2 / K_Ga)

W = feedback(C * F_1 * F_2, 1 / K_r);

figure, step(W / dcgain(W))
figure, bode(W / dcgain(W))

% Valutazioni extra
u = K_c * m_d^2 / m_i^2

% Discretizzazione
w_b = 3.97;
T_1 = (2 * pi) / (20 * w_b)
T = T_1;
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

F = (3 * (s + 2)) / (s * (s^3 + 6.5 * s^2 + 12 * s + 4.5));
N = 10;

[m_G, m_phi, w_mG, w_mphi] = margin(F);
w_pi = 1.5943;
[module, phase] = bode(F, w_pi)

% Calcolo dei parametri per il metodo Ziegler-Nichols in catena chiusa
K_P_ = m_G;
T_ = 2 * pi / w_pi;

% Calcolo dei parametri per il regolatore PID
K_P = 0.6 * K_P_
T_i = 0.5 * T_
T_d = 0.125 * T_

R = K_P * (1 + 1 / (T_i * s) + (T_d * s) / (1 + T_d / N * s));

W = feedback(R * F, 1);
figure, bode(W / dcgain(W))
