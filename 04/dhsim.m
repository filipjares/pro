
%% Inicializace robotického toolboxu
% addpath rvctools/
% startup_rvc.m

%% Definice robotu s použitím D-H parametrů

clear L;
clear r;
%            theta  d      a      alpha
L(1) = Link([0      0      0     -pi/2]);
L(2) = Link([0      0     -0.280  0]);
L(3) = Link([0      0     -0.1      pi/2]);
L(4) = Link([0      0.315  0     -pi/2]);
L(5) = Link([0      0      0      pi/2]);
L(6) = Link([0      0.085  0      0]);

% V aktualni verzi Robotics Toolboxu je prikaz 'robot' nahrazen tridou 'SerialLink'
r = SerialLink(L, 'name', 'HW-04');

% vektor nulových úhlů a vektor úhlů robotu v "pracovní" poloze
qz = [0 0 0 0 0 0];
qr = [0 pi/3 pi/6 0 pi/4 0];

% trajektorie přechodu
t = [0:.056:2]';
q = jtraj(qz, qr, t);

%% Robot v "nulové" poloze
figure(2); clf;
r.plot(qz);
axis equal tight;

%% Robot v "pracovní" poloze
figure(1); clf;
r.plot(qr);
axis equal tight;

%% Animace přechodu robotu z "nulové" polohy do "pracovní" polohy
figure(3); clf;
r.plot(qz);
axis equal tight;
r.plot(q);



