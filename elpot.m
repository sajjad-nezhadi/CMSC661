function elpot()
% Solves grad * (a(x,y) * grad u) = 0 on [0,1]^2\(B union C)
% BC's: u(0,y) = u(1,y), u_x(0,y) = u_x(1,y); 
% u_y(x,0) = u_y(x,1) = 0;
% u = 0 in B, u = 1 in C


n = 401;
xaux = linspace(0,1,n);
[x,y] = meshgrid(xaux,xaux);
u = zeros(n);
ny = n;
nx = n - 1;
hx = 1/nx;
hy = 1/(ny - 1);
aaux = 2.1 + sin(2*pi*x) + cos(3*pi*y);
a = aaux(:,1 : nx);
nxy = nx*ny;

cN = 0.5*(a + circshift(a,[1 0]))/(hy^2);
cS = 0.5*(a + circshift(a,[-1 0]))/(hy^2);
cW = 0.5*(a + circshift(a,[0 1]))/(hx^2);
cE = 0.5*(a + circshift(a,[0 -1]))/(hx^2);
cP = -(cE + cW + cN + cS);
% Conversion to column vectors
cn = cN(:); cs = cS(:); cw = cW(:); ce = cE(:); cp = cP(:);

% First set up the mathrix A without taking care of B and C
A = spdiags([circshift(cw,[-ny 0]),circshift(cn,[-1 0]),cp,...
    circshift(cs,[1 0]),circshift(ce,[ny 0])],[-ny,-1,0,1,ny],nxy,nxy);
% Take care of the Neumann BCs: u_y(x,0) = u_y(x,1) = 0
for j = 1 : nx - 1 
    A(ny*j,ny*j - 1) = A(ny*j,ny*j - 1) + A(ny*j,ny*j + 1);
    A(ny*j,ny*j + 1) = 0;
    A(ny*j + 1,ny*j + 2) = A(ny*j + 1,ny*j + 2) + A(ny*j + 1,ny*j);
    A(ny*j + 1,ny*j) = 0;
end
A(1,2) = A(1,2) + cn(1);
A(nxy,nxy - 1) = A(nxy,nxy - 1) + cs(nxy);
% Take care of the periodic BCs: u(0,y) = u(1,y), u_x(0,y) = u_x(1,y);
A1 = spdiags(circshift(ce,[-nxy + ny, 0]),-nxy + ny,nxy,nxy);
A2 = spdiags(circshift(cw,[nxy - ny, 0]),nxy - ny,nxy,nxy);
A = A + A1 + A2;

% Take care of B and C
% Regions B and C are balls centered at (x1,y1) and (x2,y2) of radii r1 and
% r2 respectively
x1 = 0.3;
y1 = 0.2;
r1 = 0.1;
x2 = 0.7;
y2 = 0.8;
r2 = 0.1;
iB = find((x - x1).^2 + (y - y1).^2 < r1^2); % indices of mesh points lying in B
iC = find((x - x2).^2 + (y - y2).^2 < r2^2); % indices of mesh points lying in C

% Set the RHS due to the BC's at C
uc = zeros(nxy,1);
uc(iC) = 1;
rhs = -A*uc;
% Remove rows and columns corresponding to iBC
iBC = [iB;iC];
inBC = [1 : nxy]';
inBC(iBC) = []; % indices of mesh points not in B and not in C
A(:,iBC) = [];
A(iBC,:) = [];
rhs(iBC) = [];

% Solution
uaux = A\rhs;
ucol = zeros(nxy,1);
ucol(inBC) = uaux;
ucol(iC) = 1;
ucol(iB) = 0;
u1 = reshape(ucol,ny,nx);
u = zeros(ny,nx + 1);
u(:,1 : nx) = u1;
u(:,nx + 1) = u1(:,1);

% Plot the conductivity, the solution, and the current;
figure(1);
clf;
hold on;
ma = max(max(aaux));
mi = min(min(aaux));
contourf(xaux,xaux,aaux,linspace(mi,ma,10));
t = linspace(0,2*pi,100);
plot(x1 + r1*cos(t),y1 + r1*sin(t),'k','Linewidth',4);
plot(x2 + r2*cos(t),y2 + r2*sin(t),'k','Linewidth',4);
colorbar;
set(gca,'DataAspectRatio',[1,1,1],'Fontsize',24);
xlabel('x','Fontsize',24);
ylabel('y','Fontsize',24);
%
figure(2);
clf;
hold on;
contourf(xaux,xaux,u,[0:0.1:1]);
t = linspace(0,2*pi,100);
plot(x1 + r1*cos(t),y1 + r1*sin(t),'w','Linewidth',2);
plot(x2 + r2*cos(t),y2 + r2*sin(t),'w','Linewidth',2);
colorbar
set(gca,'DataAspectRatio',[1,1,1],'Fontsize',24);
xlabel('x','Fontsize',24);
ylabel('y','Fontsize',24);

%
figure(3);
clf;
hold on;
curx = zeros(ny,nx + 1);
cury = zeros(ny,nx + 1);
curx(:,1 : nx) = a.*(0.5*(circshift(u1,[0 -1]) - circshift(u1,[0 1]))/hx);
cury(:,1 : nx) = a.*(0.5*(circshift(u1,[-1 0]) - circshift(u1,[1 0]))/hy);
curx(:,nx + 1) = curx(:,1);
cury(:,nx + 1) = cury(:,1);
cury(1,:) = 0;
cury(ny,:) = 0;
curx = -curx;
cury = -cury;
ecur = sqrt(curx.^2 + cury.^2);
ecmax = max(max(ecur));
contourf(xaux,xaux,ecur,linspace(0,ecmax,10));
plot(x1 + r1*cos(t),y1 + r1*sin(t),'w','Linewidth',2);
plot(x2 + r2*cos(t),y2 + r2*sin(t),'w','Linewidth',2);
ix = [1 : 8 : nx];
iy = [1 : 8 : ny];
ac = sqrt(curx(iy,ix).^2 + cury(iy,ix).^2 + 1e-12);
quiver(x(iy,ix),y(iy,ix),curx(iy,ix)./ac,cury(iy,ix)./ac,'color','w')
contour(xaux,xaux,aaux,linspace(mi,ma,10),'color','k');
colorbar
set(gca,'DataAspectRatio',[1,1,1],'Fontsize',24);
xlabel('x','Fontsize',24);
ylabel('y','Fontsize',24);
caxis([mi ma]);
end


