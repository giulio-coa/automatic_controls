close all
clear all

% Dati
s = tf('s');

F = (10 * (s + 10)) / (s^2 + 0.5 * s + 25);
d_1 = 1;
d_2 = 1;
d_3 = 1;
w_d3 = 1000;

% Analisi delle specifiche statiche
K_F = dcgain(F);

h = 1
K_c = 2000

% Analisi delle specifiche dinamiche
t_s = 0.04;
w_bdes = 3 / t_s;
w_cdes = 0.63 * w_bdes

s_max = 0.35;
M_r = (1 + s_max) / 0.85;
m_phi = 60 - 5 * 20 * log10(M_r)
m_phi = 50;

% Prima approssimazione della funzione d'anello
G_a1 = (K_c * F) / s^h;
[module_1, phase_1] = bode(G_a1, w_cdes);
recovery = m_phi - phase_1 - 180;

% Rete derivatrice
m_d = 4;
x_d = sqrt(m_d);
tau_d = x_d / w_cdes;
R_d = (1 + tau_d * s) / (1 + tau_d / m_d * s);

R_d = R_d^2;

G_a2 = G_a1 * R_d;
[module_2, phase_2] = bode(G_a2, w_cdes);

% Rete integratrice
m_i = 6.5;
x_i = 150;
tau_i = x_i / w_cdes;
R_i = (1 + tau_i / m_i * s) / (1 + tau_i * s);

R_i = R_i^2;

G_a3 = G_a2 * R_i;
[module_3, phase_3] = bode(G_a3, w_cdes);

% Controllore
C = (K_c * R_d * R_i) / s^h

G_a = C * F;
K_Ga = dcgain(s^h * G_a);

figure, margin(G_a)

% Verifica delle specifiche
e_r = abs(1 / K_Ga)

e_d2 = abs(d_2 / K_Ga)

W = feedback(C * F, 1);

figure, step(W / dcgain(W))

% Valutazioni extra
figure, bode(W / dcgain(W))

w_0 = 0.5;
amplitude = 1;
W_e = feedback(1, C * F);
[module_e, phase_e] = bode(W_e, w_0);
e_max = module_e * amplitude

W_u = feedback(C, F);
[module_u, phase_u] = bode(W_u, w_d3);
u_dp = d_3 * module_u

% Discretizzazione
w_b = 71.1;
T_1 = (2 * pi) / (20 * w_b)
T = 0.001;
G_aZOH = G_a / (1 + T / 2 * s);
C_z = c2d(C, T, 'tustin')

figure, margin(G_aZOH)

F_z = c2d(F, T, 'tustin');
W_z = feedback(C_z * F_z, 1);

figure, step(W_z / dcgain(W_z))
