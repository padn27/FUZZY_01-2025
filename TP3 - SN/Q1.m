%% MODELO TS EXATO

clc; clear; close all;

% Parâmetros físicos
p = 0.5; c = 1.0; g = 9.81; v = 5.0;

% Número de regras TS
N = 8;
sector_bounds = linspace(-pi/6, pi/6, N+1);

% Inicia células para A_i e B_i
A_TS = cell(1,N);
B_TS = cell(1,N);

for i = 1:N
    a = sector_bounds(i);
    b = sector_bounds(i+1);
    
    % Valores médios para sin e cos em cada setor
    s = (sin(a) + sin(b))/2;
    coss = (cos(a) + cos(b))/2;
    
    % Matrizes A e B do sistema linearizado no setor i
    A = zeros(4,4);
    B = zeros(4,1);
    
    A(1,2) = v;
    A(3,4) = 1;
    A(4,2) = (1)*coss * v^2 / p; % aproximação x_2 * sin(x_3) = 0
    A(4,3) = (g * s + (1)*(-s) * v^2) / p;
    
    B(2) = 1;
    B(4) = c * coss * v / p;
    
    A_TS{i} = A;
    B_TS{i} = B;
end

disp('Modelo TS exato com 8 regras criado:');
for i=1:N
    fprintf('Regra %d:\n', i);
    disp('A = '); disp(A_TS{i});
    disp('B = '); disp(B_TS{i});
end
