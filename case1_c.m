
q0 = 1000 * 0.01;  % heat source at centre of plate
k = 1;      % thermal conductivity
ne = 10;    % number of elements
nn = 11;    % number of nodes
syms e n
syms T [nn 1]

temp = sym(zeros([1 nn]));
temp(nn) = 1-e-n;
pos = zeros([nn 2]);
pos(nn, :) = [1/2 1/2];

F1 = zeros([nn 1]);
K1 = zeros(nn);

% Element-level
for i = 1:ne
    e0 = 0;
    n0 = 0;       % location of heat source in natural coordinates
    N = temp;     % shape functions
    N(i) = e;
    if i == ne
        N(1) = n;
    else
        N(i+1) = n;
    end
    
    P = pos;      % position of nodes in cartesian coordinates
    if i == 1
        P(2, :) = [1 0];
    elseif i > 1 && i < 5
        P(i, :) = [1 (i-2)/3];
        P(i+1, :) = [1 (i-1)/3];
    elseif i > 4 && i < 8
        P(i, :) = [(8-i)/3 1];
        P(i+1, :) = [(7-i)/3 1];
    elseif i > 7 && i < 10
        P(i, :) = [0 (11-i)/3];
        P(i+1, :) = [0 (10-i)/3];
    else
        P(i, :) = [0 (11-i)/3];
        P(1, :) = [0 0];
    end

    x = N * P(:, 1);
    y = N * P(:, 2);
    
    Q = [diff(N, e);
        diff(N, n)];
    
    J = [diff(x, e) diff(y, e); diff(x, n) diff(y, n)];
    iJ = J^-1;
    I1 = k * Q' * (iJ' * iJ) * Q * det(J);
    E1 = (q0/ne) * (subs(N, [e n], [e0 n0]))';
    E2 = int(int(I1, n, 0, 1-e), e, 0, 1);

    F1 = F1 + E1;
    K1 = K1 + E2;

end

% Assembly
T(1) = 298;
T(2) = 298;
fixed = [1 2];          % Nodes for which temperature is fixed (input in increasing node number only)
free = [3 4 5 6 7 8 9 10 11];

K = K1(free, free);
F = F1(free) - K1(free, fixed) * [T(1); T(2)];

T(free) = K^-1 * F


% Plotting
X = zeros([11 1]);
Y = zeros([11 1]);

X(1) = 0;
Y(1) = 0;
for i = 2:5
    X(i) = 1;
    Y(i) = (i-2)/3;
end
for i = 6:8
    X(i) = (8-i)/3;
    Y(i) = 1;
end
for i = 9:10
    X(i) = 0;
    Y(i) = (11-i)/3;
end
X(11) = 1/2;
Y(11) = 1/2;

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