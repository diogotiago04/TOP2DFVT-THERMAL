%THERMAL ANALYSIS OF 2D HOMOGENEOUS SURFACES BY FINITE-VOLUME THEORY
clear;
clc;

%__________________________________________________________USER-DEFINED
H = 1; L = 1;                                               % volume dimensions
nx2 = 50; nx1 = 50;                                         % subvolumes in x and X2
k = 25;                                                     % thermal conductivity of the material
px1 = 35;                                                   % Subvolume for the temperature field plot along the X1 direction
px2 = 1251;                                                 % Subvolume for the temperature field plot along the X2 direction
tempFI = 100;                                               % temperature of the bottom faces
tempFD = 200;                                               % temperature of the right faces
tempFS = 100;                                               % temperature of the top faces
tempFE = 100;                                               % temperature of the left faces

%_________________________________________________________INITIALIZATION OF DESIGN VARIABLE VECTOR
l = L/nx1; h = H/nx2;                                       % subvolume dimensions
FI = 1:1:nx1;                                               % bottom faces
FS = FI + nx1*nx2;                                          % top faces
FE = 1 +(nx2+1)*nx1:(nx1+1):1 +(nx2+1)*nx1+(nx1+1)*(nx2-1); % left faces
FD = FE + nx1;                                              % right faces

%_________________________________________________________AUXILIARY MATRICES EMPLOYED IN THE FINITE-VOLUME THEORY
a = ones(4,1);
N1 = [0,-1]; N2 = [1,0]; N3 = [0,1]; N4 = [-1,0];
N = [N1,zeros(1,6);zeros(1,2),N2,zeros(1,4);
    zeros(1,4),N3,zeros(1,2); zeros(1,6),N4];
A = [0 -h/2 0 h^2/4; l/2 0 l^2/4 0;
     0 h/2 0 h^2/4; -l/2 0 l^2/4 0];
E = [0 0 0 0; 0 -1 0 3*h/2; -1 0 -3*l/2 0; 0 0 0 0;
     0 0 0 0; 0 -1 0 -3*h/2; -1 0 3*l/2 0; 0 0 0 0];
k1 = [k 0 0 0 0 0 0 0; 0 k 0 0 0 0 0 0; 0 0 k 0 0 0 0 0; 0 0 0 k 0 0 0 0;
      0 0 0 0 k 0 0 0; 0 0 0 0 0 k 0 0; 0 0 0 0 0 0 k 0; 0 0 0 0 0 0 0 k];
B = N*k1*E;
ab = (B*A^-1*a\B*A^-1);
Ab = A^-1*(eye(4)-a*ab);
ke = B*Ab;
ke = [ke(1,:)*l;ke(2,:)*h;ke(3,:)*l;ke(4,:)*h];

%_________________________________________________________GLOBAL THERMAL CONDUCTIVITY MATRIX
Nx2 = nx1*(nx2+1);
Nx1 = nx2*(nx1+1);
Nf = Nx2+Nx1;
KG = sparse(Nf,Nf);
[i,j] = ndgrid(1:nx1,1:nx2);
q = (i+(j-1)*nx1);
faces = [q(:),q(:)+j(:)+nx1*(nx2+1),q(:)+nx1,q(:)+j(:)-1+nx1*(nx2+1)];
for i=1:nx1*nx2
    KG(faces(i,:),faces(i,:)) = KG(faces(i,:),faces(i,:))+ke;
end

%_________________________________________________________FACES OF THE MATERIAL AFTER DISCRETIZATION
facestemp = [FI,FD,FS,FE];                                 % faces with prescribed temperatures
temperatureFaces = zeros(1,length(facestemp));
facesfree = 1:max(faces(:));
facesfree = facesfree(~ismember(facesfree,facestemp));     % faces without prescribed temperatures                       

%_________________________________________________________PRESCRIPTION OF TEMPERATURES ON THE BOTTOM, RIGHT, TOP, AND LEFT FACES
for i=1:length(FI)
    temperatureFI(1,i) = tempFI;
end
for i=1:length(FD)
    temperatureFD(1,i) = tempFD;
end
for i=1:length(FS)
    temperatureFS(1,i) = tempFS;
end
for i=1:length(FE)
    temperatureFE(1,i) = tempFE;
end
temperatureFaces = [temperatureFI,temperatureFD,temperatureFS,temperatureFE]';

%_________________________________________________________CALCULATION OF TEMPERATURE ON FREE FACES
Qa = sparse(length(facesfree),1);
Kaa = KG(facesfree,facesfree);   
Kab = KG(facesfree,facestemp); 
Ta = (Kaa^-1)*(Qa-Kab*temperatureFaces); 
T = zeros(Nf,1);  
T(facesfree,1) = Ta; 
T(facestemp,1) = temperatureFaces; 
Tij = zeros(nx1*nx2,4);
for b=1:nx1*nx2
    sv(b).T = T(faces(b,:)); 
    T00(b) = ab*sv(b).T;
    sv(b).T00 = ab*sv(b).T; 
    sv(b).Tij = A^-1*sv(b).T - A^-1*a*sv(b).T00; 
end

%_________________________________________________________AVERAGE TEMPERATURE AT THE NODES
Tno = zeros(((nx1*nx2)*4),3);
count = 1;
for j = 1:nx2
    for i = 1:nx1 
        q = i + (j-1)*nx1;
        vert = [(i-1)*l,(j-1)*h; i*l,(j-1)*h; i*l,j*h; (i-1)*l,j*h];
        if q==px1
            cood_px = vert(1,1);
        end
        if q==px2
            cood_py = vert(1,2);
        end
        
        t = nodaltemp(sv(q),l,h);
        
               if count<=(nx2*nx1)*4
                Tno(count,1) = vert(1,1);
                Tno(count,2) = vert(1,2);
                Tno(count,3) = t(1,1);
                Tno(count+1,1) = vert(2,1);
                Tno(count+1,2) = vert(2,2);
                Tno(count+1,3) = t(1,2);
                Tno(count+2,1) = vert(3,1);
                Tno(count+2,2) = vert(3,2);
                Tno(count+2,3) = t(1,3);
                Tno(count+3,1) = vert(4,1);
                Tno(count+3,2) = vert(4,2);
                Tno(count+3,3) = t(1,4);
                count = count + 4;
               end
    end
end
count = 0;
[coords_u, ~, dx] = unique(Tno(:, 1:2), 'rows', 'stable'); % identifies the unique node coordinates in the matrix Tno
temperatures = accumarray(dx, Tno(:, 3), [], @sum);        % sums the temperatures associated with each unique coordinate
count = accumarray(dx, 1);                                 % counts the number of occurrences of each coordinate
me_temperatures = temperatures ./ count;                   % calculates the average temperature for each coordinate
Tno_M = [coords_u, me_temperatures];                       % matrix with coordinates of each node and their average temperatures

%_________________________________________________________PLOTTING OF THE NODAL TEMPERATURE FIELD

xnodes = unique(Tno_M(:,1));   
ynodes = unique(Tno_M(:,2));   
[~,~,ix] = unique(Tno_M(:,1));   
[~,~,iy] = unique(Tno_M(:,2));   
Zsum = accumarray([iy ix], Tno_M(:,3), [], @sum);
Zcnt = accumarray([iy ix], 1,          [], @sum);
Z    = Zsum ./ Zcnt;                   

[Xg, Yg] = meshgrid(xnodes, ynodes);
figure(1);
surf(Xg, Yg, Z);
view(2);                
shading flat;          
axis equal tight;
colormap(turbo);
colorbar;
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 16, ...
       'Color', 'black', 'FontWeight', 'bold');
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 16, ...
       'Color', 'black', 'FontWeight', 'bold');
title('Finite-Volume Theory (FVT)', 'Interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize',12);

   
%_________________________________________________________PLOTTING OF THE TEMPERATURE FIELD WITH A CUT ALONG THE X2 AXIS
x = 0:0.01:L;
Temp = zeros(size(x));
for dx = 1:length(x)
    i = x(dx);
    for j = 1:2:100
        denom = sinh((j*pi*L)/H);
        if denom == 0 || isinf(denom) || isnan(denom)
            k(j) = 0;
        else
            k(j) = ((400/(j*pi*denom)) * sinh(j*pi*i/H) * sin(j*pi*cood_py/H));
        end
    end
    T = 100 + sum(k);
    Temp(dx) = T;
end

ppx = l/2:l:L;
for i=1:1:length(ppx)
    Tpx(i)=sv(px2).T(1);
    px2=px2+1;
end
figure(2)
plot(x, Temp, 'k');
hold on;
plot(ppx,Tpx,'ro');
hold off;
legend({'Analytical Solution', 'FVT'}, 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black', 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black', 'FontWeight', 'bold')
title1 = sprintf('Temperature field at $x_2 = %.2f$',cood_py);
title(title1, 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black', 'FontWeight', 'bold');
ylim([floor(min(Tpx))-10 floor(max(Tpx))+10]);
yticks(100:20:floor(max(Tpx))+5);
xticks(0:0.2:L);
xlim([0 L]);

%_________________________________________________________PLOTTING OF THE TEMPERATURE FIELD WITH A CUT ALONG THE X1 AXIS
x = 0:0.01:H;
Temp = zeros(size(x));
for dx = 1:length(x)
    i = x(dx);
    for j = 1:2:100
        denom = sinh((j*pi*L)/H);
        if denom == 0 || isinf(denom) || isnan(denom)
            k(j) = 0;
        else
            k(j) = ((400/(j*pi*denom)) * sinh(j*pi*cood_px/H) * sin(j*pi*i/H));
        end
    end
    T = 100 + sum(k);
    Temp(dx) = T;
end

ppy = h/2:h:H;
for i=1:1:length(ppy)
    Tpy(i)=sv(px1).T(4);
    px1=px1+nx1;
end
figure(3)
plot(x, Temp, 'k');
hold on;
plot(ppy,Tpy,'ro');
hold off;
legend({'Analytical Solution', 'FVT'}, 'Interpreter', 'latex');
xlabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black', 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black', 'FontWeight', 'bold')
title1 = sprintf('Temperature field at $x_1 = %.2f$',cood_px);
title(title1, 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black', 'FontWeight', 'bold');
ylim([floor(min(Tpy))-10 floor(max(Tpy))+10]);
yticks(100:20:floor(max(Tpy))+5);
xticks(0:0.2:H);
xlim([0 H]);

%_________________________________________________________FUNCTION FOR CALCULATING NODAL TEMPERATURE
function T=nodaltemp(SV,l,h)
Tij = SV.Tij;
T00 = SV.T00;

vert = [-l/2,-h/2; l/2,-h/2; l/2,h/2; -l/2,h/2];

for i=1:4 
    x1 = vert(i,1);
    x2 = vert(i,2);
   T(i) = T00 + x1*Tij(1) + x2*Tij(2) + (0.5*(3*(x1^2) - ((l^2)/4))*Tij(3))...
       + (0.5*(3*(x2^2) - ((h^2)/4))*Tij(4));
end
end
