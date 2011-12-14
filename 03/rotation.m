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
axis equal tight;

%% 5.
% Nakreslete bod Y do kterého se pohne bod X.

Y = R'*X + o - o;

plot_vector(o, Y, 'k', 'y''');
hold off;
grid on;

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

% (1) osu pohybu (R,o=[0;0;0]) znázorníme jako 'r';
% (2) osu pohybu (R,o=[1;1;1]) jako 'a'
% Protoze maji pohyby (1) a (2) stejnou matici R, budou jejich osy pohybu
% rovnobezne

% Pohyb (1)
I = eye(3);
r = null(R' - I);

hold on;
plot_vector([0 0 0]', r, 'r', 'r');

% Pohyb (2) -- pro (R,o=[1;1;1]):
o_beta = -o;
a = null((R' - I)^2);                   % == lambda * null(R' - I)
p = pinv((R' - I)^2)*(R' - I)*o_beta;   % bod na ose pohybu
plot_vector(p - 2*a, 4*a, 'g', 'a');

%% 3. Najděte osu rotace r a nakreslete ji. Jaký má osa rotace vztah k ose pohybu a?

% Osa rotace pohybu (1) i (2) prochazi bodem O = [0 0 0]' a ma smerovy
% vektor odpovidajici r. (V pripade pohybu (1) je osa rotacetotozna s osou
% pohybu)

r = null(R' - I);

% Osa rotace je rovnobezna k ose pohybu
    
%% 4. Pozorovaný vztah dokažte. (Použijte rovnici (R'-I)^2*x_\beta = (R'-I)*o_\beta.)

% Osa rotace:
% pro poby na ose rotace x, plati:
%     R'*x == x,                neboli
%     (R' - I)*x == [0 0 0]'
%
% Vektory na ose rotace (a tedy i smerovy vektor osy rotace -- osa rotace
% prochazi nulou) patri do jadra matice (R' - I)
%
% Osa pohybu:
% Vektory xi na ose pohybu splnuji rovnici
%     ((R' - I)^2)*xi == (R' - I)*o_beta.
% Vezmeme-li rozdil dvou ruznych bodu lezicich
% na ose pohybu x1, x2, mame z predchoziho vztahu:
%     ((R' - I)^2)*(x1 - x2) == (R' - I)*o_beta - (R' - I)*o_beta == 0.
% Smerove vektory osy pohybu tedy patri do jadra matice (R' - I)^2
%
% Pro rotacni matici R a jednotkovou matici I plati, ze v pruseciku
% jadra a oboru hodnot matice (R' - I) lezi pouze nulovy vektor.
% Z toho vyplyva, ze jadro matice (R' - I) je stejne jako jadro
% matice (R' - I)^2. Z toho vyplyva, ze smerove vektory osy rotace
% a osy pohybu lezi ve stejnem jednorozmernem podprostoru a osa
% pohybu a osa rotace jsou tedy rovnobezne.


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
