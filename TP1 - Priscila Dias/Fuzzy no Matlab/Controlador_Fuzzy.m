%% Limpeza do ambiente
clear; clc; close all;

%% Parâmetros da simulação
dt = 0.1;
tempo = 0:dt:60;

% Trajetória circular do líder
raio = 5;
velocidade_angular_lider = 0.1;
x_l = raio * cos(velocidade_angular_lider * tempo);
y_l = raio * sin(velocidade_angular_lider * tempo);

% Estado inicial do robô seguidor
x_s = zeros(size(tempo)); y_s = zeros(size(tempo)); theta_s = zeros(size(tempo));
x_s(1) = 0; y_s(1) = -8; theta_s(1) = 0;

%% Criação do controlador fuzzy
fis = mamfis('Name','RoboSeguidor');

% Entradas
fis = addInput(fis, [-pi pi], 'Name', 'ErroAngular');
fis = addMF(fis, 'ErroAngular', 'trimf', [-pi -pi -pi/2], 'Name', 'NG');
fis = addMF(fis, 'ErroAngular', 'trimf', [-pi -pi/2 0], 'Name', 'NP');
fis = addMF(fis, 'ErroAngular', 'trimf', [-1 0 1], 'Name', 'Z');
fis = addMF(fis, 'ErroAngular', 'trimf', [0 pi/2 pi], 'Name', 'PP');
fis = addMF(fis, 'ErroAngular', 'trimf', [pi/2 pi pi], 'Name', 'PG');

fis = addInput(fis, [0 10], 'Name', 'ErroPosicao');
fis = addMF(fis, 'ErroPosicao', 'trimf', [0 0 2], 'Name', 'Perto');
fis = addMF(fis, 'ErroPosicao', 'trimf', [1.5 3.5 6], 'Name', 'Medio');
fis = addMF(fis, 'ErroPosicao', 'trimf', [5 7.5 10], 'Name', 'Longe');

% Saídas
fis = addOutput(fis, [-2 2], 'Name', 'VelAng');
fis = addMF(fis, 'VelAng', 'trimf', [-1.5 -1.5 0], 'Name', 'Negativa');
fis = addMF(fis, 'VelAng', 'trimf', [-0.5 0 0.5], 'Name', 'Zero');
fis = addMF(fis, 'VelAng', 'trimf', [0 1.5 1.5], 'Name', 'Positiva');

fis = addOutput(fis, [0 5], 'Name', 'VelLin');
fis = addMF(fis, 'VelLin', 'trimf', [0 0 1], 'Name', 'Lenta');
fis = addMF(fis, 'VelLin', 'trimf', [0.5 1 1.5], 'Name', 'Media');
fis = addMF(fis, 'VelLin', 'trimf', [1 2 2], 'Name', 'Rapida');

% Regras
ruleList = [
    1 1 1 1 1 1;
    2 1 1 1 1 1;
    3 1 2 2 1 1;
    4 1 3 2 1 1;
    5 1 3 1 1 1;
    1 2 1 2 1 1;
    2 2 1 2 1 1;
    3 2 2 3 1 1;
    4 2 3 2 1 1;
    5 2 3 2 1 1;
    1 3 1 3 1 1;
    2 3 1 3 1 1;
    3 3 2 3 1 1;
    4 3 3 3 1 1;
    5 3 3 3 1 1;
];
fis = addRule(fis, ruleList);

%% Simulação do seguidor
for t = 1:length(tempo)-1
    xd = x_l(t); yd = y_l(t);
    ep = sqrt((xd - x_s(t))^2 + (yd - y_s(t))^2);
    theta_d = atan2(yd - y_s(t), xd - x_s(t));
    ea = wrapToPi(theta_d - theta_s(t));

    entrada = [ea ep];
    out = evalfis(fis, entrada);
    v = out(2); w = out(1);

    x_s(t+1) = x_s(t) + dt * v * cos(theta_s(t));
    y_s(t+1) = y_s(t) + dt * v * sin(theta_s(t));
    theta_s(t+1) = theta_s(t) + dt * w;
end

%% Plotagem das trajetórias
figure;
plot(x_l, y_l, 'b--', 'LineWidth', 2); hold on;
plot(x_s, y_s, 'r-', 'LineWidth', 2);
legend('Líder (Circular)', 'Seguidor (Fuzzy)');
xlabel('x [m]'); ylabel('y [m]');
title('Trajetória do Robô Líder e Seguidor');
grid on; axis equal;

%% Visualização do Controlador Fuzzy
writeFIS(fis, 'controlador_seguidor.fis');
fuzzyLogicDesigner('controlador_seguidor.fis');

figure;
plotfis(fis);

figure;
gensurf(fis, [1 2], 1); % VelAng
figure;
gensurf(fis, [1 2], 2); % VelLin
