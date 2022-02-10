function hw1()
n = 100;
nx = round(2*pi*n);
ny = 2*n;
tx = linspace(0,2*pi,nx);
ty = linspace(0,2,ny);
[x,y] = meshgrid(tx,ty);

f = ((x > pi/2) & (x < 3*pi/2)) .* -cos(x);

f1 = f(1 : ny - 1,1 : nx - 1);
f_aux = f1(:);
h = 2/(nx - 1);
u = zeros(ny,nx);

% Set up the matrix A
I = speye(ny - 1);
Ix = speye(nx - 1);
e = ones(nx-1,1);
T = spdiags([e -4*e e],[-1:1],nx-1,nx-1);
T(1,nx-1) = 1;
T(nx-1,1) = 1;
S = kron(spdiags([e e],[-1 1],ny-1,ny-1),Ix);
l = zeros(ny-1);
l(1,2) = 1;
P = kron(l,Ix);
A = (kron(I,T) + S + P)/h^2;
    
% Solve the linear system
u_aux = A\f_aux;
u(1:ny-1,1:nx-1) = reshape(u_aux,ny-1,nx-1);

% plot the solution
figure(1);
clf; hold on; grid;
ma = max(max(u));
mi = min(min(u));
contourf(x,y,u,linspace(mi,ma,20));
title('numerical solution');
xlabel('x');
ylabel('y');

end
