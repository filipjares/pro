%% Část odpovídající druhému domácímu úkolu, níže následuje třetí úkol
% 1. Simulujte pohyb se zadanou maticí R a translací o, kde o má
% význam o_{\beta'}, jako v rovnici 3.8 v [PRO-2011-Lecture-02.pdf].
% Uvažujte bázi \beta jako standardní bázi.

% matice rotace, priblizna
R = [0.8047   -0.5059   -0.3106; ...
     0.3106    0.8047   -0.5059; ...
     0.5059    0.3106    0.8047];

% matice rotace, lepe:
[U,D,V] = svd(R);
R = U*V';

% translace;
o = [1;1;1]; % zadane mam o = o_{\beta'}
o = -R'*o;   % chci aby o = o'_{\beta}

%% 2.
% Najděte souřadnice vektorů \beta' v \beta a naopak.
% Nakreslete vektory \beta a \beta' ve standardní bázi.

% Odpověďi:
% - souřadnice vektorů \beta v \beta' jsou sloupce matice R
% - souřadnice vektorů \beta' v \beta jsou sloupce transponované matice R
%   (resp. řádky matice R)

% pocatek a souradne soustav beta a beta' (beta2 = beta');
beta = [ ...
    1 0 0; ...
    0 1 0; ...
    0 0 1];
beta2 = R';

figure(1); % nakreslíme souřadné soustavy ve standardní bázi:
plot_coordinate_frame([0 0 0]', beta, 'k', 'beta', 'O'); hold on;
plot_coordinate_frame([0 0 0]', beta2, 'b', 'beta'''); hold off;
grid on;
axis equal;

%% 3.Nakreslete souřadné soustavy (O=0,\beta) a (O',\beta').

figure(2);
plot_coordinate_frame([0 0 0]', beta, 'k', 'beta', 'O'); hold on;
plot_coordinate_frame(o, beta2, 'b', 'beta''', 'O''');
grid on;
axis equal;

%% 4.
% Nakreslete polohové vektory bodu X = [1;2;3]
% vzhledem k (O=0,\beta) a (O',\beta').

X = [1 2 3]';
X2 = X - o;

plot_vector([0 0 0]', X, 'g', 'x');
plot_vector(o, X2, 'm', 'x''');

%% 5.
% Nakreslete bod Y do kterého se pohne bod X.

Y = R'*X + o - o;

plot_vector(o, Y, 'k', 'y''');
hold off;
grid on;
axis equal;

%% Část odpovídající třetímu domácímu úkolu následuje:
%
%    1. Uvažujte stejný pohyb jako v HW-02 a pokračujte v kreslení do obrázku z HW-02.
%    2. Najděte osy pohybů (R,o=[0;0;0]) a (R,o=[1;1;1]). Nakreslete je do obrázku.
%    3. Najděte osu rotace r a nakreslete ji. Jaký má osa rotace vztah k ose pohybu a?
%    4. Pozorovaný vztah dokažte. (Použijte rovnici (R'-I)^2*x_\beta = (R'-I)*o_\beta.)
%    5. Najděte rovinu p kolmou na osu rotace, množinu jejích generátorů (vektory, které
%       ji generují), a nakreslete generátory.
%    6. Jaký je vztah mezi generátory p a maticí (R'-I)?
%    7. Najděte bod P, ve kterém osa pohybu a protíná rovinu p, a nakreslete jej.
%    8. Nakreslete bod P', který vznikne otočením bodu P, a bod P'', který vznikne
%       posunutím bodu P'. Jaký je vztah P, P' a osy a?
%    9. Jaký je vztah mezi osou rotace r a osou pohybu a když:
%          1. R' = I
%          2. o_\beta = 0
%          3. o_\beta je vlastním vektorem R'

%% 2. Najděte osy pohybů (R,o=[0;0;0]) a (R,o=[1;1;1]). Nakreslete je do obrázku.
% osu pohybu (R,o=[0;0;0]) znázorníme jako 'r'; osu pohybu (R,o=[1;1;1]) jako 'a'

I = eye(3);
r = null(R' - I);

hold on;
plot_vector([0 0 0]', r, 'r', 'r');

% pro (R,o=[1;1;1]):
o_beta = -o;
a = null((R' - I)^2);
p = pinv((R' - I)^2)*(R' - I)*o_beta;
plot_vector(p - 2*a, 4*a, 'g', 'a');

%% 3. Najděte osu rotace r a nakreslete ji. Jaký má osa rotace vztah k ose pohybu a?

r = null(R' - I);

% osa rotace je rovnoběžná k ose pohybu
    
%% 4. Pozorovaný vztah dokažte. (Použijte rovnici (R'-I)^2*x_\beta = (R'-I)*o_\beta.)
%% 5. Najděte rovinu p kolmou na osu rotace, množinu jejích generátorů (vektory, které
%     ji generují), a nakreslete generátory.

p = orth((R' - I));
plot_vector([0 0 0]', p(:,1), 'b', 'p_1');
plot_vector([0 0 0]', p(:,2), 'b', 'p_2');
hold off;
%% 6. Jaký je vztah mezi generátory p a maticí (R'-I)?

% generátory p, které generují rovinu p kolmou na osu rotace odpovídají
% bázi oboru hodnot matice (R'-I)

%% 7. Najděte bod P, ve kterém osa pohybu a protíná rovinu p, a nakreslete jej.
%% 8. Nakreslete bod P', který vznikne otočením bodu P, a bod P'', který vznikne
%     posunutím bodu P'. Jaký je vztah P, P' a osy a?
%% 9. Jaký je vztah mezi osou rotace r a osou pohybu a když:
%       1. R' = I
%       2. o_\beta = 0
%       3. o_\beta je vlastním vektorem R'
