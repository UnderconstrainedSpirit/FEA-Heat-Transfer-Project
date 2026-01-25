q0 = 1000 * 0.01;  % heat source at centre of plate
k = 1;      % thermal conductivity
ne = 100;    % number of elements
nn = 101;    % number of nodes
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
    if i > 0 && i < 26
        P(i, :) = [(i-1)/25 0];
        P(i+1, :) = [i/25 0];
    elseif i > 25 && i < 51
        P(i, :) = [1 (i-26)/25];
        P(i+1, :) = [1 (i-25)/25];
    elseif i > 50 && i < 76
        P(i, :) = [(76-i)/25 1];
        P(i+1, :) = [(75-i)/25 1];
    elseif i > 75 && i < 100
        P(i, :) = [0 (101-i)/25];
        P(i+1, :) = [0 (100-i)/25];
    else
        P(i, :) = [0 (101-i)/25];
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

fixed = zeros([1 26]);     % Nodes for which temperature is fixed (input in increasing node number only)
free = zeros([1 75]);
known = zeros([26 1]);
for iter = 1:26
    fixed(iter) = iter;
    T(iter) = 298;
    known(iter) = 298;
end

for iter2 = 1:75
    free(iter2) = iter2 + 26;
end

K = K1(free, free);
F = F1(free) - K1(free, fixed) * known;

T(free) = K^-1 * F


% Plotting

X = zeros([101 1]);
Y = zeros([101 1]);
for i = 1:26
    X(i) = (i-1)/25;
    Y(i) = 0;
end
for i = 27:51
    X(i) = 1;
    Y(i) = (i-26)/25;
end
for i = 52:76
    X(i) = (76-i)/25;
    Y(i) = 1;
end
for i = 77:100
    X(i) = 0;
    Y(i) = (101-i)/25;
end
X(101) = 1/2;
Y(101) = 1/2;

Z = double(T);

nq = 50;  % number of grid points in each direction
[xq,yq] = meshgrid(linspace(0,1,nq), linspace(0,1,nq));

Zq = griddata(X, Y, Z, xq, yq);

% draw filled contours
figure;
contourf(xq, yq, Zq, 20, 'LineColor','none');
colorbar;
axis square tight;
xlabel('x'); ylabel('y');
title('Temperature Contours');


x_line = linspace(0,1,200);
y_line = 0.5*ones(size(x_line));

T_line = griddata(X, Y, Z, x_line, y_line);

figure;
plot(x_line, T_line, 'LineWidth',1.5);
xlabel('x at y=0.5');
ylabel('T(x,0.5)');
grid on;
title('Temperature Variation Along Centerline');
