q0 = 1000 * 0.01;  % heat source at centre of plate
k = 1;      % thermal conductivity
syms e n
syms T [5 1]

% Element 1
e0 = 0;
n0 = 0;       % location of heat source in natural coordinates
N = [e n 0 0 1-e-n];    % shape functions
P = [0 0; 1 0; 0 0; 0 0; 1/2 1/2];      % positions of nodes in cartesian coordinates
x = N * P(:, 1);
y = N * P(:, 2);

Q = [diff(N, e);
    diff(N, n)];

J = [diff(x, e) diff(y, e); diff(x, n) diff(y, n)];
iJ = J^-1;
I1 = k * Q' * (iJ' * iJ) * Q * det(J);
E11 = (q0/4) * (subs(N, [e n], [e0 n0]))';
E12 = int(int(I1, n, 0, 1-e), e, 0, 1);

% Element 2
e0 = 0;
n0 = 0;       % location of heat source in natural coordinates
N = [0 e n 0 1-e-n];    % shape functions
P = [0 0; 1 0; 1 1; 0 0; 1/2 1/2];      % positions of nodes in cartesian coordinates
x = N * P(:, 1);
y = N * P(:, 2);

Q = [diff(N, e);
    diff(N, n)];

J = [diff(x, e) diff(y, e); diff(x, n) diff(y, n)];
iJ = J^-1;
I2 = k * Q' * (iJ' * iJ) * Q * det(J);
E21 = (q0/4) * (subs(N, [e n], [e0 n0]))';
E22 = int(int(I2, n, 0, 1-e), e, 0, 1);
 
%Element 3
e0 = 0;
n0 = 0;       % location of heat source in natural coordinates
N = [0 0 e n 1-e-n];    % shape functions
P = [0 0; 0 0; 1 1; 0 1; 1/2 1/2];      % positions of nodes in cartesian coordinates
x = N * P(:, 1);
y = N * P(:, 2);

Q = [diff(N, e);
    diff(N, n)];

J = [diff(x, e) diff(y, e); diff(x, n) diff(y, n)];
iJ = J^-1;
I3 = k * Q' * (iJ' * iJ) * Q * det(J);
E31 = (q0/4) * (subs(N, [e n], [e0 n0]))';
E32 = int(int(I3, n, 0, 1-e), e, 0, 1);

%Element 4
e0 = 0;
n0 = 0;       % location of heat source in natural coordinates
N = [n 0 0 e 1-e-n];    % shape functions
P = [0 0; 0 0; 0 0; 0 1; 1/2 1/2];      % positions of nodes in cartesian coordinates
x = N * P(:, 1);
y = N * P(:, 2);

Q = [diff(N, e);
    diff(N, n)];

J = [diff(x, e) diff(y, e); diff(x, n) diff(y, n)];
iJ = J^-1;
I4 = k * Q' * (iJ' * iJ) * Q * det(J);
E41 = (q0/4) * (subs(N, [e n], [e0 n0]))';
E42 = int(int(I4, n, 0, 1-e), e, 0, 1);

% Assembly
T(1) = 298;
T(2) = 298;
fixed = [1 2];     % Nodes for which temperature is fixed (input in increasing node number only)
free = [3 4 5];

K1 = E12 + E22 + E32 + E42;
F1 = (E11 + E21 + E31 + E41);

K = K1(free, free);
F = F1(free) - K1(free, fixed) * [T(1); T(2)];

T(free) = K^-1 * F

% Plotting

X = [0 1 1 0 1/2]';
Y = [0 0 1 1 1/2]';
Z = double(T);

nq = 50;
[xq, yq] = meshgrid(linspace(0, 1, nq), linspace(0, 1, nq));
Zq = griddata(X, Y, Z, xq, yq);

figure
contourf(xq, yq, Zq, 20, 'LineColor', 'none');
colorbar;
axis square tight;
xlabel('x'); ylabel('y');
title("Temperature Contour");


x_line = linspace(0, 1, 200);
y_line = 0.5 * ones(size(x_line));

T_line = griddata(X, Y, Z, x_line, y_line, 'Linear');

figure
plot(x_line, T_line);
xlabel('x');
ylabel('T(x, 0.5)');
grid on;
title("Temperature variations along centerline")