
%% Inicializace robotického toolboxu
% addpath rvctools/
% startup_rvc

%% Definice D-H parametru robotu

a1 = 85;
a2 = 280;
a3 = 100;
a4 = 0;
a5 = 0;
a6 = 0;
d1 = 350;
d2 = 0;
d3 = 0;
d4 = 315;
d5 = 0;
d6 = 85;
mu1 = -1;
mu2 = 0;
mu3 = -1;
mu4 = 1;
mu5 = -1;
mu6 = 0;
lambda1 = 0;
lambda2 = 1;
lambda3 = 0;
lambda4 = 0;
lambda5 = 0;
lambda6 = 1;

alpha1 = -pi/2;
alpha2 = 0;
alpha3 = -pi/2;
alpha4 = pi/2;
alpha5 = -pi/2;
alpha6 = 0;

%% Definice robotu s pouzitim D-H parametru

clear L;
clear r;
%            theta  d      a      alpha
L(1) = Link([0      d1	   a1     alpha1]);
L(2) = Link([0      d2	   a2     alpha2]);
L(3) = Link([0      d3	   a3     alpha3]);
L(4) = Link([0      d4	   a4     alpha4]);
L(5) = Link([0      d5	   a5     alpha5]);
L(6) = Link([0      d6	   a6     alpha6]);

% V aktualni verzi Robotics Toolboxu je prikaz 'robot' nahrazen tridou 'SerialLink'
r = SerialLink(L, 'name', 'HW-07 #03');
r1 = SerialLink(L, 'name', 'HW-07 #03');

%%

load('P0.mat', '-ascii');   % pozice odpovidajici kloubovym souradnicim t0
load('P1.mat', '-ascii');   % pozice odpovidajici kloubovym souradnicim t1
load('P2.mat', '-ascii');   % pozice odpovidajici kloubovym souradnicim t2
load('P3.mat', '-ascii');   % pozice odpovidajici kloubovym souradnicim t3 == t0
load('t0.mat', '-ascii');   % kloubove souradnice t0 == [0 0 0 0 0 0]'
load('t1.mat', '-ascii');   % kloubove souradnice t1
load('t2.mat', '-ascii');   % kloubove souradnice t2
load('t3.mat', '-ascii');   % kloubove souradnice t3 == [0 0 0 0 0 0]'

% Transformacni matice ziskane z predchoziho volanim
% ve tvaru          MhV0 = RobCoordsToMhVMatrix(P0);
load('MhV0.mat', '-ascii');
load('MhV1.mat', '-ascii');
load('MhV2.mat', '-ascii');


%%

    figure(1); % nakreslíme souřadné soustavy ve standardní bázi:
    plot_coordinate_frame([0 0 0]', 500*eye(3), 'k', 'x', 'world'); hold on;
    plot_coordinate_frame(P0(1:3,1), 500*anglesToMtx(P0(4:6,1)), 'b', 'x', 'eff_0');
    plot_coordinate_frame(P1(1:3,1), 500*anglesToMtx(P1(4:6,1)), 'r', 'x', 'eff_1');
    plot_coordinate_frame(P2(1:3,1), 500*anglesToMtx(P2(4:6,1)), 'm', 'x', 'eff_2');
    hold off;
    grid on;
    axis equal;

%% Reseni IKU

load('sol0.mat', '-ascii');
load('sol1.mat', '-ascii');
load('sol2.mat', '-ascii');

%%

for i = 1:size(sol0, 2)
    q = sol0(:,i);

    figure(i);
    r.plot(q');
    axis equal tight;
    view([-1 -0.5 1]);
end

q1 = sol1(:,1);
figure(5);
r1.plot(q1');
axis equal tight;
view([-1 -0.5 1]);

