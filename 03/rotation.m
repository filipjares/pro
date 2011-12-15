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
plot3(X(1), X(2), X(3), 'ob', 'markerfacecolor', 'b');
plot_vector(o, X2, 'm', 'x''');
axis equal tight;

%% 5.
% Nakreslete bod Y do kterého se pohne bod X.

Y = R'*X + o;

plot_vector(o, Y - o, 'k', 'y''');
plot3(Y(1), Y(2), Y(3), 'ob', 'markerfacecolor', 'b');
hold off;
grid on;
axis tight;

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
% Protože mají pohyby (1) a (2) stejnou matici R, budou jejich osy pohybu
% rovnoběžné

% Pohyb (1)
I = eye(3);
r = null(R' - I);

hold on;
plot_vector([0 0 0]', r, 'r', 'r');

% Pohyb (2) -- pro (R,o=[1;1;1]):
o_beta = -o;
a = null((R' - I)^2);                   % == lambda * null(R' - I)
P = pinv((R' - I)^2)*(R' - I)*o_beta;   % bod na ose pohybu
plot_vector(P - 3*a, 4*a, 'g', 'a');

%% 3. Najděte osu rotace r a nakreslete ji. Jaký má osa rotace vztah k ose pohybu a?

% Osa rotace pohybu (1) i (2) prochází bodem O = [0 0 0]' a má směrový
% vektor odpovídající r. (V případě pohybu (1) je osa rotace totožná s osou
% pohybu)

r = null(R' - I);

% Osa rotace je rovnoběžná k ose pohybu
    
%% 4. Pozorovaný vztah dokažte. (Použijte rovnici (R'-I)^2*x_\beta = (R'-I)*o_\beta.)

% Osa rotace:
% pro body na ose rotace x, platí:
%     R'*x == x,                neboli
%     (R' - I)*x == [0 0 0]'
%
% Vektory na ose rotace (a tedy i směrový vektor osy rotace -- osa rotace
% prochází nulou) patří do jádra matice (R' - I)
%
% Osa pohybu:
% Vektory xi na ose pohybu splňují rovnici
%     ((R' - I)^2)*xi == (R' - I)*o_beta.
% Vezmeme-li rozdíl dvou různých bodů ležících
% na ose pohybu x1, x2, máme z předchozího vztahu:
%     ((R' - I)^2)*(x1 - x2) == (R' - I)*o_beta - (R' - I)*o_beta == 0.
% Směrové vektory osy pohybu tedy patří do jádra matice (R' - I)^2
%
% Pro rotační matici R a jednotkovou matici I platí, že v průsečíku
% jádra a oboru hodnot matice (R' - I) leží pouze nulový vektor.
% Z toho vyplývá, že jádro matice (R' - I) je stejné jako jádro
% matice (R' - I)^2. A tedy směrové vektory osy rotace
% a osy pohybu leží ve stejném jednorozměrném podprostoru a osa
% pohybu a osa rotace jsou tedy rovnoběžné.


%% 5. Najděte rovinu p kolmou na osu rotace, množinu jejích generátorů (vektory, které
%     ji generují), a nakreslete generátory.

% Generátory roviny p jsou vektory p_1 a p_2. Rovina p je tvořena všemi
% jejich lineárními kombinacemi.

p = orth((R' - I));     % baze oboru hodnot matice (R' - I)
plot_vector([0 0 0]', p(:,1), 'b', 'p_1');
plot_vector([0 0 0]', p(:,2), 'b', 'p_2');

%% 6. Jaký je vztah mezi generátory p a maticí (R'-I)?

% generátory p, které generují rovinu p kolmou na osu rotace odpovídají
% bázi oboru hodnot matice (R'-I)

%% 7. Najděte bod P, ve kterém osa pohybu a protíná rovinu p, a nakreslete jej.

% bod P jsme uz nalezli v ramci bodu 2.
plot3(P(1), P(2), P(3), 'ob', 'markerfacecolor', 'b');
text(P(1), P(2), P(3), 'P');


%% 8. Nakreslete bod P', který vznikne otočením bodu P, a bod P'', který vznikne
%     posunutím bodu P'. Jaký je vztah P, P' a osy a?

P2 = R'*P;          % P'
plot3(P2(1), P2(2), P2(3), 'ob', 'markerfacecolor', 'b');
text(P2(1), P2(2), P2(3), 'P''');

P3 = P2 - o_beta;   % P''
plot3(P3(1), P3(2), P3(3), 'ob', 'markerfacecolor', 'b');
text(P3(1), P3(2), P3(3), 'P''''');

% Spojnice bodů P a P' je kolmá na osu a. Posunutí bodu P' vrátí bod P
% do polohy P'', t.j. na osu pohybu. Pohyb přesouvá body na ose pohybu
% na body na ose pohybu (P na P'').

hold off;

%% 9. Jaký je vztah mezi osou rotace r a osou pohybu a když:
%       1. R' = I
%       2. o_\beta = 0
%       3. o_\beta je vlastním vektorem R'

% 1) Jedná se o pohyb, který je tvořen pouze posunutím. Za osu rotace
% bychom v tomto případě mohli považovat libovolnou přímku procházející
% počátkem. Libovolnou přímku v prostoru můžeme považovat za osu pohybu.
%
% 2) Pohyb tvořený pouze rotací - osa rotace je totožná s osou pohybu a
% body na ose pohybu zůstávají na svém místě.
%
% 3) V tomto případě je pohyb tvořen rotací a posunutím rovnoběžným
% s osou rotace a osa pohybu je opět totožná s osou pohybu, v tomto případě
% ovšem body na ose pohybu nezůstávají na svém místě, ale posouvají se ve
% směru o\beta.
