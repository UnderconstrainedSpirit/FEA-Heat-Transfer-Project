q0 = 1000 * 0.01;  % heat source at centre of plate
k = 1;      % thermal conductivity
syms e n
syms T [4 1]

% Element 1
e0 = 0;
n0 = 1/2;       % location of heat source in natural coordinates
N1 = e;
N2 = n;
N3 = 0;
N4 = 1-e-n;
x = n;
y = 1-e-n;

Q = [diff(N1, e) diff(N2, e) diff(N3, e) diff(N4, e); diff(N1, n) diff(N2, n) diff(N3, n) diff(N4, n)];
J = [diff(x, e) diff(y, e); diff(x, n) diff(y, n)];
iJ = J^-1;
I1 = k * Q' * (iJ' * iJ) * Q * det(J);
E11 = (q0) * [e0; n0; 0; 1-e0-n0];
E12 = int(int(I1, n, 0, 1-e), e, 0, 1);

% Element 2
e0 = 1/2;
n0 = 0;
N1 = 0;
N2 = e;
N3 = n;
N4 = 1-e-n;
x = e+n;
y = 1-e;

Q = [diff(N1, e) diff(N2, e) diff(N3, e) diff(N4, e); diff(N1, n) diff(N2, n) diff(N3, n) diff(N4, n)];
J = [diff(x, e) diff(y, e); diff(x, n) diff(y, n)];
iJ = J^-1;
I2 = k * Q' * (iJ' * iJ) * Q * det(J);
E21 = (q0) * [0; e0; n0; 1-e0-n0];
E22 = int(int(I2, n, 0, 1-e), e, 0, 1);

% Assembly
given = [1, 2];     % Nodes for which temperature is fixed (input in increasing node number only)
T(1) = 298;
T(2) = 298;
fixed = [1 2];
free = [3 4];

K1 = E12 + E22;
F1 = (E11 + E21);

K = K1(free, free);
F = F1(free) - K1(free, fixed) * [T(1); T(2)];

T(free) = K^-1 * F

% Plotting

X = [0 1 1 0]';
Y = [0 0 1 1]';
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

