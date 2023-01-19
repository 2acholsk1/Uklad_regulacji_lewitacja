%%
clc;
clear all;
close all;

%% parametry obiektu i stale
nr_albumu = [1 4 7 5 8 9];

m = 0.023/(1+nr_albumu(5)+nr_albumu(6)); % kg
g = 9.81; % m/s^2
FemP1 = 1.7521*10^-2; % H
FemP2 = 5.8231*10^-3; % m
f1 = 1.4142*10^-4; % m*s
f2 = 4.5626*10^-3; % m
ki = 2.5165; % A
ci = 0.0243*sign(nr_albumu(2)-nr_albumu(6)); % A
d = 0.066*(1 + 0.1*nr_albumu(6)); % m
bd = 0.06; % m
xd = d-bd; % m

%% punkt pracy
x10 = 0.008; % m
x20 = 0; % m/s
x30 = 0.2; % A % <--- TU W POLECENIU JEST 0.75, dla 0.2 JEST OK
u10 = 1/ki*(x30-ci); % [-]

%% obl a_mn
a21 = 1/2/m * x30^2 * FemP1/FemP2^2 * exp(-x10/FemP2);
a23 = -1/m * x30 * FemP1/FemP2 * exp(-x10/FemP2);
a31 = 1/f1 * exp(x10/f2) * (ki*u10 + ci - x30);
a33 = -f2/f1 * exp(x10/f2);
b31 = ki * f2/f1 * exp(x10/f2);

%% ctrb/obsv/stab
A = [0 1 0; a21 0 a23; a31 0 a33];
B = [0 0; 0 g; b31 0];
C = [1 0 0];
D = 0;
sys = ss(A,B,C,D);

disp("sterowalnosc; obserwowalnosc; stabilnosc");
rank(ctrb(sys)) == rank(A)
rank(obsv(sys)) == rank(A)
isstable(sys)

%% wyznaczanie transmitancji
G = tf(sys);
G = G(1)
Gf = 1/G;

%% nastawy PID
kp = 110;
ki = 30;
kd = 4;

%% sprzezenie od stanu
% u = x_des - Kx
% x' = Ax + B(x_des - Kx)
% x' = Ax + Bx_des - BKx
% x' = (A - BK)x + Bx_des

K = place(A, B(:, 1), [-28, -60, -90])
A_cl = A - B(:, 1)*K
B_cl = B

%% obserwator stanu
Tp = 0.00001;
p = [-600, -700, -600];
L = acker(A', C', p)';

%% postac diagonalna
csys = canon(sys)

%% postac sterowalna
syms s
sI = s*eye(3);
Ms = det(sI-A);
pretty(Ms);

As = [0 1 0; 0 0 1; 3.8138*10^5 2.0473*10^3 -186.2891]
Bs = [0;0;1]

S = [B(:,1) A*B(:,1) A^2*B(:,1)];
Ss = [Bs As*Bs As^2*Bs];

P = Ss*S^-1;

Cs = C*P^-1
Ds = D
