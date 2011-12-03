
EPS = 100*eps;

% silenost odpovidajici volani randn('seed', 0) ve starsich verzich matlabu
s = RandStream('mcg16807','Seed', 0)
RandStream.setDefaultStream(s)

% 1. nahodne vygenerovana matice A
% randn('seed', #);
A = randn(3);
assert(rank(A) == 3);

% 2. Vypoctete svd rozklad matice A
[U, D, V] = svd(A);
assert(norm(U*D*V' - A) < EPS);

% 3. Najdete matici Abar, hod(Abar) == 2, aby norm(Abar - A) byla minimalni
%    a vypoctete norm(Abar-A) pomoci SVD
Dbar = D;
Dbar(3,3) = 0;
Abar = U*Dbar*V';
assert(rank(Abar) == 2);

[~,De,~] = svd(Abar-A);
assert(De(1,1) - norm(Abar-A) < EPS);

% 4. Perturbujte Abar a zkonstruujte B = Abar + tau*randn(3), pricemz
%    tau nabyva hodnot 1e-3, 1e-6, 1e-9

tau = [1e-3, 1e-6, 1e-9];
for i = 1:size(tau,2)
    E{i} = tau(i)*randn(3);
    B{i} = Abar + E{i};
    assert(rank(B{i}) == 3);
end;

% 5. Najdete bazi reseni soustav Abar*x = 0, rank(Abar) == 2 a 
%    Bbar*x = 0, hod(Bbar) == 2

% baze reseni Abar*x = 0 je vektor
vb = V(:,3);
assert((vb'*vb - 1) < EPS);

for i = 1:size(tau,2)
    % najdeme Bbar
    [Ub{i}, Db{i}, Vb{i}] = svd(B{i});
    Dbbar{i} = Db{i}; Dbbar{i}(3,3) = 0;
    Bbar{i} = Ub{i}*Dbbar{i}*Vb{i}';
    assert(rank(Bbar{i}) == 2);
    % baze reseni Bbar*x = 0 (ma dimenzi 1):
    [Ui, Di, Vi] = svd(Bbar{i});
    % baze reseni:
    vB{i} = Vi(:,3);
    assert((vB{i}'*vb - 1) < EPS);
    % rozdil mezi bazi vb{i} a vb{i} (vektory muzou ovsem byt opacne)
    vBdiff{i} = min(norm(vB{i} - vb), norm(vB{i} - (-vb)));
end
    
% 6. Najdete vektor na kterem se realizuje norma Abar, tzn vektor x,
%    norm(x) == 1, pro nejz je norm(Abar) == norm(Abar*x)

% vektor, na kterem se realizuje norma Abar je vektor v1:
[~, ~, VV] = svd(Abar);
v1 = VV(:,1);
% ukazeme, ze je to pravda:
assert(norm(Abar*v1) - norm(Abar) < EPS); 

% 7. Porovnejte E a (B - Bbar) a jejich normy

for i = 1:size(tau,2)
   disp(i);
   disp('E, Bbar a jejich normy');
   E{i}
   B{i} - Bbar{i}
   norm(E{i})
   norm(B{i}-Bbar{i})  % je (musi!) byt mensi nez norma E
end

% 8. C = randn(6,3) -> SVD, pozorujte velikot matic

C = randn(6,3);
[Uc,Dc,Vc] = svd(C);


% 9. Vypoctete inv(A) z U,D,V bez použití inv

% 10. b = Abar*randn(3,1). Vyreste pomoci SVD soustavy A*x = b, Abar*x = b

% 11. x na kruznici -> nakreslete x, A*x, Abar*x, B*x, *Bbar*x v R^2