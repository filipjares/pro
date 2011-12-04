
% kontrolovat budeme pomoci assertu srovnavajicich s EPS:
EPS = 100*eps;

% silenost odpovidajici volani randn('seed', 0) ve starsich verzich matlabu
s = RandStream('mcg16807','Seed', 0);
RandStream.setDefaultStream(s);

% 1. nahodne vygenerovana matice A
% randn('seed', #);
A = randn(3);
assert(rank(A) == 3);
disp('1.: Nahodne vygenerovana matice A s hodnosti 3:');
disp(A);

% 2. Vypoctete svd rozklad matice A
[U, D, V] = svd(A);
assert(norm(U*D*V' - A) < EPS);

% 3. Najdete matici Abar, hod(Abar) == 2, aby norm(Abar - A) byla minimalni
%    a vypoctete norm(Abar-A) pomoci SVD
Dbar = D;
Dbar(3,3) = 0;
Abar = U*Dbar*V';
assert(rank(Abar) == 2);
disp('3.: Matice Abar, sestavena z matice A tak, aby hod(Abar) == 2) a aby norm(Abar-A) byla minimalni:');
disp(Abar);

[~,De,~] = svd(Abar-A);
assert(De(1,1) - norm(Abar-A) < EPS);
disp('norma matice Abar odpovida prvku D(1,1) matice D z SVD rozkladu matice Abar (U*D*V'' = Abar):');
disp(De(1,1));

% 4. Perturbujte Abar a zkonstruujte B = Abar + tau*randn(3), pricemz
%    tau nabyva hodnot 1e-3, 1e-6, 1e-9

tau = [1e-3, 1e-6, 1e-9];
for i = 1:size(tau,2)
    E{i} = tau(i)*randn(3);
    B{i} = Abar + E{i};
    assert(rank(B{i}) == 3);
end;
disp('4.: K prvkum matice Abar byly pridany nahodne poruchy a tak vznikly matice B{i}, i = 1,2,3');
disp(' '); % existuje lepsi zpusob, jak odradkovat?

% 5. Najdete bazi reseni soustav Abar*x = 0, rank(Abar) == 2 a 
%    Bbar*x = 0, hod(Bbar) == 2

vb = V(:,3);
assert((vb'*vb - 1) < EPS);
disp('5.: Baze reseni Abar*x = 0 je vektor V(:,3) z SVD rozkladu matice Abar = U*D*V'':');
disp(vb);

for i = 1:size(tau,2)
    % najdeme Bbar tak, ze nahradime prvek (3,3) v diagonalni matici
    % z SVD rozkladu matice B nulovym prvkem
    [Ub{i}, Db{i}, Vb{i}] = svd(B{i});
    Dbbar{i} = Db{i}; Dbbar{i}(3,3) = 0;
    Bbar{i} = Ub{i}*Dbbar{i}*Vb{i}';
    assert(rank(Bbar{i}) == 2);
    % baze reseni Bbar*x = 0 (ma dimenzi 1):
    [Ui, Di, Vi] = svd(Bbar{i});
    % (jednoprvkove) baze vB{i} nulovych prostoru jednotlivych matic Bbar{i}:
    vB{i} = Vi(:,3);
    assert((vB{i}'*vb - 1) < EPS);
    % rozdil mezi bazi vb{i} a vb{i} (vektory muzou ovsem byt opacne)
    vBdiff{i} = min(norm(vB{i} - vb), norm(vB{i} - (-vb)));
end
    
% 6. Najdete vektor na kterem se realizuje norma Abar, tzn vektor x,
%    norm(x) == 1, pro nejz je norm(Abar) == norm(Abar*x)

[~, ~, VV] = svd(Abar);
v1 = VV(:,1);
% ukazeme, ze je to pravda:
assert(norm(Abar*v1) - norm(Abar) < EPS);
disp('6.: Vektor, na kterem se realizuje norma Abar je takovy jednotkovy vektor v,');
disp('    pro nejz je norm(''Abar*v'') nejvetsi. Takovy vektor v odpovida prvnimu');
disp('    sloupci matice V z SVD rozkladu matice Abar = U*D*V'':');
disp(v1);

% 7. Porovnejte E a (B - Bbar) a jejich normy

disp('7.: "Porovnejte E a (B - Bbar) a jejich normy:');
for i = 1:size(tau,2)
   %disp('E, Bbar a jejich normy');
   % E{i};
   % B{i} - Bbar{i};
   % B{i} - Bbar{i} je (musi byt!) mensi nez norma E
   disp([ 'norm(E' num2str(i) ') = ' num2str(norm(E{i})) '; norm(B' num2str(i) ' - Bbar' num2str(i) ') = ' num2str(norm(B{i}-Bbar{i})) ]);
end
disp(' '); % existuje lepsi zpusob, jak odradkovat?

% 8. C = randn(6,3) -> SVD, pozorujte velikot matic

C = randn(6,3);
[UC, DC, VC] = svd(C);
disp('8.: C = randn(6,3) -> SVD, pozorujte velikot matic. C == U*D*V'':');
disp('    size(C) = ');
disp(size(C))
disp('    size(U) = ');
disp(size(UC));
disp('    size(D) = ');
disp(size(DC));
disp('    size(V) = ');
disp(size(VC));

% 9. Vypoctete inv(A) z U,D,V bez použití inv

disp('9.: Vypoctete inv(A) z U,D,V bez použití inv.'); disp(' ');
disp('SVD rozklad matice A je: A == U*D*V''. Muzeme proto psat:');
disp('inv(A) == inv(V'')*inv(D)*inv(U) == (V'')''*iD*U'' == V*iD*U'',');
disp('kde iD je matice stejnych rozmeru jako matice D, ale s prvky');
disp('odpovidajicimi prevracenym hodnotam prvku v D. Tedy inv(A) = ');
iD = diag(diag(1./D));
iA = V*iD*U';
disp(iA);
assert(norm(inv(A) - iA) < EPS);

% 10. b = Abar*randn(3,1). Vyreste pomoci SVD soustavy A*x = b, Abar*x = b

disp('10.: b = Abar*randn(3,1). Vyreste pomoci SVD soustavy A*x = b, Abar*x = b');
disp(' '); % existuje lepsi zpusob, jak odradkovat?

b = Abar*randn(3,1);
disp('Rovnice A*x = b je ekvivalentni rovnici A*x - b = 0, takze muzeme vyuzit reseni soustavy');
disp('C*y = 0, kde y = [x; 1] a C = [A -b]. K reseni teto soustavy pozijeme SVD stejne jako');
disp('v bode 5. pri hledani baze nuloveho prostoru');
disp(' '); % existuje lepsi zpusob, jak odradkovat?
disp('Pripad pro rovnici A*x = b: C = [A -b] = ');
C = [A -b];
disp(C);
[UC1, DC1, VC1] = svd(C);
% baze nuloveho prostoru matice C ma dimenzi 1 (protoze rank(C) == rank(A) == 3 a 4 - 3 = 1)
% reseni soustavy A*x = b je jedine
assert(norm(C*VC1(:,end) - zeros(3,1)) < EPS);
assert(rank(C) == rank(A) & rank(A) == 3);
y = VC1(:,end);
x = VC1(1:3,end)/VC1(4,end);
disp('Resenim soustavy A*x = b je jediny vektor x = ');
disp(x);

disp('Pripad pro rovnici Abar*x = b: C = [Abar -b] = ');
C2 = [Abar -b];
disp(C2);
[UC2, DC2, VC2] = svd(C2);
% baze nuloveho prostoru matice C2 ma dimenzi 2 (protoze rank(C) == rank(A) == 2 a 4 - 2 = 2)
assert(norm(C2*VC2(:,end) - zeros(3,1)) < EPS);
assert(norm(C2*VC2(:,end-1) - zeros(3,1)) < EPS);
assert(rank(C2) == rank(Abar) & rank(Abar) == 2);
y1 = VC2(:,end-1);
y2 = VC2(:,end);
yy = y1 - y1(4)/y2(4)*y2;
yy = yy(1:3);
y1 = y1(1:3)/y1(4);
y2 = y2(1:3)/y2(4);
disp('Resenim soustavy Abar*x = b je kazdy vektor ve tvaru x = x0 + k*x1 pro libovolne k.');
disp('pritom, x0 je partikularni reseni a je napr. x0 =');
disp(y2);
disp('a x1 je reseni pridruzene homogenni soustavy Abar*x1 = 0 a je x1 = ');
disp(yy);
for i=1:100
    k = randn(1);
    assert(norm(Abar*(y2 + k*yy) - b) < 100*EPS);
end

% 11. Pro body x na kruznici nakreslete spolu se samotnymi x jejich obrazy
%     A*x, Abar*x, B*x, *Bbar*x v R^2

% zkonstruujeme potrebne matice:
A2 = randn(2);
assert(rank(A2) == 2);
[U2, D2, V2] = svd(A2);
D2bar = D2; D2bar(2,2) = 0;
A2bar = U2*D2bar*V2';
for i = 1:size(tau,2)       % tau pouzijeme stejne jako vyse
    E2{i} = tau(i)*randn(2);
    B2{i} = A2bar + E2{i};
    assert(rank(B2{i}) == 2);
    % najdeme Bbar tak, ze nahradime prvek (2,2) v diagonalni matici
    % z SVD rozkladu matice B nulovym prvkem
    [U2b{i}, D2b{i}, V2b{i}] = svd(B2{i});
    D2bbar{i} = D2b{i}; D2bbar{i}(2,2) = 0;
    B2bar{i} = U2b{i}*D2bbar{i}*V2b{i}';
    assert(rank(B2bar{i}) == 1);
end;

% a vsechno nakreslime
figure(1);
t = [0:0.1:2*pi];
t2 = [0.05:0.1:(2*pi+0.05)];
X = cos(t); X2 = cos(t2);
Y = sin(t); Y2 = sin(t2);
plot(X, Y, 'r.');               % x
axis equal;
hold on;
XX = A2*[X;Y];
plot(XX(1,:), XX(2,:), 'y.');   % A*x
XX = A2bar*[X2;Y2];
plot(XX(1,:), XX(2,:), 'bx');   % Abar*x
XX = B2{1}*[X;Y];
plot(XX(1,:), XX(2,:), 'mo', 'MarkerSize', 8);   % B2*x (staci nakreslit pro jedine tau)
XX = B2bar{1}*[X;Y];
plot(XX(1,:), XX(2,:), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'g');   % B2bar*x
hold off;
% a k obrazku jeste popis
xlabel('x');
ylabel('y');
legend('x', 'A*x', 'Abar*x', 'B2*x', 'B2bar*x');
title('Znazorneni bodu x, A*x, Abar*x, B*x, *Bbar*x v rovine');
