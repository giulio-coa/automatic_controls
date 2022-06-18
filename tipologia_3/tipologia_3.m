close all
clear all

% Dati
s = tf('s');

A = 9;
G_p = -0.65 / (s * (s^2 + 4 * s + 1.75));
T_p = 1;
A_1 = 5.5 * 10^-3;
A_2 = 5.5 * 10^-3;
A_p = 10^-3;
w_p = 30;

% Analisi delle specifiche statiche
K_Gp = dcgain(s * G_p);

h = 0
K_c = -1.5

% Analisi delle specifiche dinamiche
t_s = 1;
w_bdes = 3 / t_s;
w_cdes = 0.63 * w_bdes
w_b = w_bdes;
w_c = 2;

s_max = 0.3;
M_r = (1 + s_max) / 0.85;
M_rdB = 20 * log10(M_r);
m_phi = 60 - 5 * M_rdB
m_phi = 42;

% Prima approssimazione della funzione d'anello
G_a1 = (K_c * A * G_p) / (s^h * T_p);
[module_1, phase_1] = bode(G_a1, w_c);
recovery = m_phi - phase_1 - 180;

% Rete derivatrice
m_d = 14;
x_d = sqrt(m_d);
tau_d = x_d / w_c;
R_d = (1 + tau_d * s) / (1 + tau_d / m_d * s);

G_a2 = G_a1 * R_d;
[module_2, phase_2] = bode(G_a2, w_c);

% Rete integratrice
m_i = 2;
x_i = 80;
tau_i = x_i / w_c;
R_i = (1 + tau_i / m_i * s) / (1 + tau_i * s);

G_a3 = G_a2 * R_i;
[module_3, phase_3] = bode(G_a3, w_c);

% Controllore
C = (K_c * R_d * R_i) / s^h

G_a = 1 / T_p * C * A * G_p;
K_Ga = dcgain(s^(h + 1) * G_a);

figure, margin(G_a)

% Verifica delle specifiche
e_r = abs(1 / K_Ga)

e_d1 = abs(A_1 / dcgain(C * A))

e_d2 = abs(A_2 / K_Ga)

W = feedback(C * A * G_p, T_p);

figure, step(W / dcgain(W))

% Valutazioni extra
figure, bode(W / dcgain(W))

W_u = feedback(C, A * G_p * T_p);
[module_u, phase_u] = bode(W_u, w_p);
u_dp = A_p * module_u

% Discretizzazione
w_b = 3.5;
T_1 = (2 * pi) / (20 * w_b)
T = 0.02;
G_aZOH = G_a / (1 + T / 2 * s);
C_z = c2d(C, T, 'tustin')

figure, margin(G_aZOH)

G_pz = c2d(G_p, T, 'tustin');
W_z = feedback(C_z * A * G_pz, T_p);

figure, step(W_z / dcgain(W_z))
