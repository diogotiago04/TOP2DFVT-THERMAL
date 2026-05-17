    %THERMAL ANALYSIS OF 2D HOMOGENEOUS SURFACES BY FINITE-VOLUME THEORY
clear;
clc;
%__________________________________________________________USER-DEFINED
H = 1; L = 1;                                               % volume dimensions
ny = 50; nx = 50;                                           % subvolumes in X1 and X2
k = 25;                                                     % thermal conductivity of the material
fx1 = 40;                                                   % Face for the temperature field plot along the X2 direction
fx2 = 1951;                                                 % Face for the temperature field plot along the X1 direction
tempFB = 100;                                               % temperature of the bottom faces
tempFR = 200;                                               % temperature of the right faces
tempFT = 100;                                               % temperature of the top faces
tempFL = 100;                                               % temperature of the left faces

%_________________________________________________________INITIALIZATION OF DESIGN VARIABLE VECTOR
l = L/nx; h = H/ny;                                         % subvolume dimensions
FB = 1:1:nx;                                                % bottom faces
FT = FB + nx*ny;                                            % top faces
FL = 1 +(ny+1)*nx:(nx+1):1 +(ny+1)*nx+(nx+1)*(ny-1);        % left faces
FR = FL + nx;                                               % right faces

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
Nx2 = nx*(ny+1);
Nx1 = ny*(nx+1);
Nf = Nx2+Nx1;
KG = sparse(Nf,Nf);
[i,j] = ndgrid(1:nx,1:ny);
q = (i+(j-1)*nx);
faces = [q(:),q(:)+j(:)+nx*(ny+1),q(:)+nx,q(:)+j(:)-1+nx*(ny+1)];
for i=1:nx*ny
    KG(faces(i,:),faces(i,:)) = KG(faces(i,:),faces(i,:))+ke;
end

%_________________________________________________________FACES OF THE MATERIAL AFTER DISCRETIZATION
facestemp = [FB,FR,FT,FL];                                 % faces with prescribed temperatures
temperatureFaces = zeros(1,length(facestemp));
facesfree = 1:max(faces(:));
facesfree = facesfree(~ismember(facesfree,facestemp));     % faces without prescribed temperatures                       

%_________________________________________________________PRESCRIPTION OF TEMPERATURES ON THE BOTTOM, RIGHT, TOP, AND LEFT FACES
for i=1:length(FB)
    temperatureFI(1,i) = tempFB;
end
for i=1:length(FR)
    temperatureFD(1,i) = tempFR;
end
for i=1:length(FT)
    temperatureFS(1,i) = tempFT;
end
for i=1:length(FL)
    temperatureFE(1,i) = tempFL;
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
Tij = zeros(nx*ny,4);
for b=1:nx*ny
    sv(b).T = T(faces(b,:)); 
    T00(b) = ab*sv(b).T;
    sv(b).T00 = ab*sv(b).T; 
    sv(b).Tij = A^-1*sv(b).T - A^-1*a*sv(b).T00; 
end

%_________________________________________________________AVERAGE TEMPERATURE AT THE NODES
Tno = zeros(((nx*ny)*4),3);
count = 1;
for j = 1:ny
    for i = 1:nx 
        q = i + (j-1)*nx;
        vert = [(i-1)*l,(j-1)*h; i*l,(j-1)*h; i*l,j*h; (i-1)*l,j*h];
        if q==fx1
            cood_px = vert(1,1);
        end
        if q==fx2
            cood_py = vert(1,2);
        end
        
        t = nodaltemp(sv(q),l,h);
        
               if count<=(ny*nx)*4
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
tt = [];
figure (1);
hold on

for j = 1:ny
    for i = 1:nx 
        q = i + (j-1)*nx;
        vert = [(i-1)*l,(j-1)*h; i*l,(j-1)*h; i*l,j*h; (i-1)*l,j*h];

        for ii = 1:4
            for jj = 1:size(Tno_M)
                if vert(ii,1) == Tno_M(jj,1) && vert(ii,2) == Tno_M(jj,2);
                    tt(ii) = Tno_M(jj,3);
                end
            end
        end

    mdl = scatteredInterpolant(vert(:,1), vert(:,2), tt', 'natural');   
    xg = linspace(min(vert(:,1)'), max(vert(:,1)'), 5);
    yg = linspace(min(vert(:,2)'), max(vert(:,2)'), 5);   
    [Xg, Yg] = meshgrid(xg, yg);
    Zg = mdl(Xg, Yg);
    surf(Xg, Yg, Zg,'edgecolor','none','facecolor','interp');    
    end
end
    colormap(turbo);
    colorbar;
    set(gca, 'FontSize', 12);
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold');
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold');
    title('Finite-Volume Theory (FVT)', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold');
   
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
    Tpx(i)=sv(fx2).T(1);
    fx2=fx2+1;
end
figure(2)
plot(x, Temp, 'k');
hold on;
plot(ppx,Tpx,'ro');
hold off;
set(gca, 'FontSize', 12);
legend({'Analytical Solution', 'FVT'}, 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold')
title1 = sprintf('Temperature field at $x_2 = %.2f$',cood_py);
title(title1, 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold');
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
    Tpy(i)=sv(fx1).T(4);
    fx1=fx1+nx;
end
figure(3)
plot(x, Temp, 'k');
hold on;
plot(ppy,Tpy,'ro');
hold off;
set(gca, 'FontSize', 12);
legend({'Analytical Solution', 'FVT'}, 'Interpreter', 'latex');
xlabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold')
title1 = sprintf('Temperature field at $x_1 = %.2f$',cood_px);
title(title1, 'Interpreter', 'latex', 'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold');
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
