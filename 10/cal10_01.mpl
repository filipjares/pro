# 
# Packages & settings
# 
restart:
with(ListTools):
with(LinearAlgebra):
with(PolynomialTools):
with(combinat, choose):
with(Groebner):
with(MatrixPolynomialAlgebra):
interface(rtablesize=24):
interface(warnlevel=0):
Digits:=30:
eps:=1e-6:
# 
# Definice procedur
# Two one-parametric motions transformatin in DH-convention(phi, theta, a,d) indexed by i
# c = cos(phi), s = sin(phi), lambda = cos(alpha), mu = sin(alpha)
dhTs := proc(i)
local M1, M2;
   M1:=Matrix(4,4,[[+cat(`c`,i),-cat(`s`,i),0,         0],
                  [ +cat(`s`,i),+cat(`c`,i),0,         0],
                  [           0,          0,1,cat(`d`,i)],
                  [           0,          0,0,         1]]);
   M2:=Matrix(4,4,[[1,        0,          0,cat(`a`,i)],  
                  [ 0,+cat(`lambda`,i),-cat(`mu`,i),         0],
                  [ 0,+cat(`mu`,i),+cat(`lambda`,i),         0],
                  [ 0,          0,          0,         1]]);
   [M1,M2];
end proc:
#
# Inverse of the DH-convention for one-aprametric DH rigid motion transformations
dhInvs := proc(M)
   local M1, M2;
   M1 := M[1];
   M2 := M[2];
   [simplify(MatrixInverse(M2),{M2[3,2]^2+M2[3,3]^2=1}),
    simplify(MatrixInverse(M1),{M1[1,1]^2+M1[2,1]^2=1})];
end proc:
#
# Rigid motion transformatin in DH-convention(phi, theta, a,d) indexed by i
# c = cos(phi), s = sin(phi), P = cos(alpha), R = sin(alpha)
dhT := proc(i)
local M;
   M:=dhTs(i);
   M[1].M[2];
end proc:
#
# Inverse of the DH-convention rigid motion transformation
dhInv := proc(M)
   simplify(MatrixInverse(M),{M[1,1]^2+M[2,1]^2=1,M[3,2]^2+M[3,3]^2=1});
end proc:
#
# Simplify using trigonometric indentities c^2+s^2=1 & lambda^2+mu^2=1
dhSimpl := proc(M,i)
   simplify(M,{cat(`c`,i)^2+cat(`s`,i)^2=1,cat(`lambda`,i)^2+cat(`mu`,i)^2=1});
end proc:
#
## Direct Kinematic Task
#
dhDKT := proc(p)
   subs(p,dhT(1).dhT(2).dhT(3).dhT(4).dhT(5).dhT(6));   
end proc:
#
# Simplify using Rotation matrin in Mh
MhSimpl := proc(M)
 simplify(
  simplify(
   simplify(
    simplify(M,
            {lx^2+ly^2+lz^2=1,mx^2+my^2+mz^2=1,nx^2+ny^2+nz^2=1}),
           {lx*mx+ly*my+lz*mz=0,lx*nx+ly*ny+lz*nz=0,mx*nx+my*ny+mz*nz=0}),
          {lx^2+mx^2+nx^2=1,ly^2+my^2+ny^2=1,lz^2+mz^2+nz^2=1}),
         {lx*ly+mx*my+nx*ny=0,lx*lz+mx*mz+nx*nz=0,lz*ly+mz*my+nz*ny=0});
end proc:
#
# Simplify a general motion matrix using rotation matrix identities in columns
rcSimp := proc(M,R)
       simplify(
        simplify(
         simplify(
          simplify(
           simplify(
            simplify(M,{R[1,1]*R[1,1]+R[2,1]*R[2,1]+R[3,1]*R[3,1]=1}),
                {R[1,1]*R[1,2]+R[2,1]*R[2,2]+R[3,1]*R[3,2]=0}),
               {R[1,1]*R[1,3]+R[2,1]*R[2,3]+R[3,1]*R[3,3]=0}),
              {R[1,2]*R[1,2]+R[2,2]*R[2,2]+R[3,2]*R[3,2]=1}),
             {R[1,2]*R[1,3]+R[2,2]*R[2,3]+R[3,2]*R[3,3]=0}),
            {R[1,3]*R[1,3]+R[2,3]*R[2,3]+R[3,3]*R[3,3]=1});
end proc:
#
# Simplify a general motion matrix using rotation matrix identities in rows
rrSimp := proc(M,R)
       simplify(
        simplify(
         simplify(
          simplify(
           simplify(
            simplify(M,{R[1,1]*R[1,1]+R[1,2]*R[1,2]+R[1,3]*R[1,3]=1}),
                       {R[1,1]*R[2,1]+R[1,2]*R[2,2]+R[1,3]*R[2,3]=0}),
                       {R[1,1]*R[3,1]+R[1,2]*R[3,2]+R[1,3]*R[3,3]=0}),
                       {R[2,1]*R[2,1]+R[2,2]*R[2,2]+R[2,3]*R[2,3]=1}),
                       {R[2,1]*R[3,1]+R[2,2]*R[3,2]+R[2,3]*R[3,3]=0}),
                       {R[3,1]*R[3,1]+R[3,2]*R[3,2]+R[3,3]*R[3,3]=1});
end proc:
#
# Matrix representation of a set of polynomials 
PolyCoeffMatrix:=proc(S,m,Ord::{ShortTermOrder, TermOrder}) 
local A,v,i,j,k,c,q;
        A:=Matrix(nops(S),nops(m),storage=sparse);
        v:=indets(m);
        for i from 1 to nops(S) do
                c:=[coeffs(expand(S[i]),v,'q')];
                q:=[q];
        	      for j from 1 to nops(m) do
                        for k from 1 to nops(q) do
                                if (m[j]=q[k]) then A[i,j]:=c[k] end if
                        end do
         	      end do
        end do;
        Matrix(A);
end proc:
#
## Cartesian product of a two lists
#
LxL:=proc(X::list,Y::list)
     Flatten(map(x->(map(y->Flatten([x,y]),Y)),X),1);
end proc:
#
## n x 1 matrix to a list conversion
#
M2L:=proc(M) 
	convert(convert(M,Vector),list);
end proc:
#
## Highlit non-zero entries
#
spy:=proc(A)
   map(x->`if`(simplify(x)=0,0, `if`(simplify(x)=1,1,`*`)) ,A):
end proc:
#
# Monomials of a set of polynomial in all indeterminates
#
PolyMonomials:=proc(S::list(ratpoly),Ord::{ShortTermOrder, TermOrder}) # Monomials of a set of polynomials
local v,m,i,c,q;
        v:=indets(S);
        m:=[];
        for i from 1 to nops(S) do
                c:=[coeffs(expand(S[i]),v,'q')];
                m:=[op(m),q];
        end do;
        m:=MakeUnique(m);
        sort(m,(t1,t2)->testorder(t2,t1,Ord));
end proc:
#
## Monomias of a set of polynomials in given indeterminates
#
PolyVarsMonomials:=proc(S::list(ratpoly),Ord::{ShortTermOrder, TermOrder}) # Monomials of a set of polynomials in variavbles of Ord
local v,m,i,c,q;
        v:={op(Ord)};
        m:=[];
        for i from 1 to nops(S) do
                c:=[coeffs(expand(S[i]),v,'q')];
                m:=[op(m),q];
        end do;
        m:=MakeUnique(m);
        sort(m,(t1,t2)->testorder(t2,t1,Ord));
end proc:
# 
# Definice procedury
IK6 := proc(Mechanism, MhV)
    #local M31,M32,M41,M42,M51,M52,iM22,iM21,iM12,iM11,Mh,iM62,iM61,
    #local ee1,ee2,l2,p2,l1,p1,pxl1,pxl2,m1x
    
    M31 :=dhTs(3)[1]:
    M32 :=dhTs(3)[2]:
    M41 :=dhTs(4)[1]:
    M42 :=dhTs(4)[2]:
    M51 :=dhTs(5)[1]:
    M52 :=dhTs(5)[2]:
    iM22:=dhInvs(dhTs(2))[1]:
    iM21:=dhInvs(dhTs(2))[2]:
    iM12:=dhInvs(dhTs(1))[1]:
    iM11:=dhInvs(dhTs(1))[2]:
    Mh := Matrix(4,4,[[lx,mx,nx,rx],[ly,my,ny,ry],[lz,mz,nz,rz],[0,0,0,1]]):
    iM62:=dhInvs(dhTs(6))[1]:
    iM61:=dhInvs(dhTs(6))[2]:
    ee1:=dhSimpl((<<1,0,0,0>|<0,-1,0,0>|<0,0,1,0>|<0,0,d2,1>>.
                    dhInv(iM22).M31.M32.M41.M42.M51.M52)[1..3,3..4],2):
    ee2:=dhSimpl((<<1,0,0,0>|<0,-1,0,0>|<0,0,1,0>|<0,0,d2,1>>.
                    iM21.iM12.iM11.Mh.iM62.iM61)[1..3,3..4],2):
    l2:=ee1[1..3,1..1]:
    p2:=ee1[1..3,2..2]:
    l1:=ee2[1..3,1..1]:
    p1:=ee2[1..3,2..2]:
    pp1:=MhSimpl(dhSimpl(dhSimpl(dhSimpl(Transpose(p1).p1,2),1),6)):
    pp2:=dhSimpl(dhSimpl(dhSimpl(dhSimpl(Transpose(p2).p2,2),3),4),5):
    pl1:=MhSimpl(dhSimpl(dhSimpl(dhSimpl(Transpose(p1).l1,2),1),6)):
    pl2:=dhSimpl(dhSimpl(dhSimpl(dhSimpl(Transpose(p2).l2,2),3),4),5):
    pxl1:=map(x->expand(x),convert(CrossProduct(convert(p1,Vector),convert(l1,Vector)),Matrix)):
    pxl2:=map(x->expand(x),convert(CrossProduct(convert(p2,Vector),convert(l2,Vector)),Matrix)):
    m1x:=map(x->expand(x),MhSimpl(dhSimpl(dhSimpl(dhSimpl(pxl1,2),1),6))):
    m2x:=map(x->expand(x),dhSimpl(dhSimpl(dhSimpl(dhSimpl(pxl2,2),3),4),5)):
    plpl1:=map(x->expand(x),ScalarMultiply(l1,pp1[1,1]) - ScalarMultiply(p1,2*pl1[1,1])):
    plpl2:=map(x->expand(x),ScalarMultiply(l2,pp2[1,1]) - ScalarMultiply(p2,2*pl2[1,1])):
    mp1:=MhSimpl(dhSimpl(dhSimpl(dhSimpl(simplify(plpl1),2),1),6)):
    mp2:=dhSimpl(dhSimpl(dhSimpl(dhSimpl(dhSimpl(simplify(plpl2),2),3),4),5),1):
    E1:=<p1,l1,pp1,pl1,m1x,mp1>:
    E2:=<p2,l2,pp2,pl2,m2x,mp2>:

    t1:=<<s1*s2,s1*c2,c1*s2,c1*c2,s1,c1,s2,c2,1>>:
    t2:=<<s4*s5,s4*c5,c4*s5,c4*c5,s4,c4,s5,c5,1>>:
    M1:=PolyCoeffMatrix(M2L(E1),M2L(t1),plex(op(indets(tE1)))):
    M2:=PolyCoeffMatrix(M2L(E2),M2L(t2),plex(op(indets(tE2)))):

    P:= <M2[1..14,1..8]|M2[1..14,9]-M1[1..14,9]>:
    Q:=  M1[1..14,1..8]:
    p:= t2:
    q:= <t1[1..8,1]>:

    MPh := {op(MhV), op(Mechanism)}:
    Qx := simplify(evalf(subs(MPh, Q))):
    Px := simplify(evalf(subs(MPh, P))):
    QxR := convert(Qx, rational, exact):
    PxR := convert(Px, rational, exact):

    U, S, Vt := SingularValues(QxR, output=[':-U', ':-S', ':-Vt']):
    A  := Transpose(U).PxR:
    A8 := A[1..8,1..9]:
    Z  := A[9..14,1..9]:

    x:=subs([s4=(2*x4)/(1+x4^2),c4=(1-x4^2)/(1+x4^2),s5=(2*x5)/(1+x5^2),c5=(1-x5^2)/(1+x5^2)],p):
    y:=ScalarMultiply(x, (1+x4^2)*(1+x5^2)):
    qy := simplify(Z.y):
    mMy := PolyVarsMonomials(M2L(qy),plex(x4,x5)):
    My := PolyCoeffMatrix(M2L(qy), mMy, plex(x4,x5)):

    qyy:=simplify(<simplify(x4*qy),qy>):
    myy:=PolyVarsMonomials(M2L(qyy),plex(x4,x5));
    Myy:=PolyCoeffMatrix(M2L(qyy),myy,plex(x4,x5)):
    Mx3:=simplify((1+x3^2)*subs([s3=(2*x3)/(1+x3^2),c3=(1-x3^2)/(1+x3^2)],Myy)):
    
    dx3:=sort(Determinant(evalf(Mx3)));
    nf:=LeadingCoefficient(dx3,tdeg(x3))^(1/12);
    dx3nf:=sort(Determinant(evalf(Mx3)/nf));
    dx3n := sort(dx3nf/LeadingCoefficient(dx3nf, tdeg(x3)));
    Cx3 := CompanionMatrix(dx3n, x3):
    x3s:=Eigenvalues(Cx3);    # solutions

    x3sR:=<select(proc(x) abs(Im(x))=0 end proc, M2L(x3s))>:
    x3s:=map(x->convert(Re(x),rational),x3sR);evalf(x3s): # real solutions

    solutions := [];
    for ixs from 1 to Dimension(x3s) do
        Mx3s:=subs(x3=x3s[ixs],Mx3):
        sv,V:= SingularValues(Mx3s, output=[':-S',':-Vt']):
        myys := Transpose(V[12..12,1..12])/V[12,12]:

        x4s:=myys[9,1]/myys[12,1]:
        x5s:=myys[11,1]/myys[12,1]:
        s3s := subs(x3=x3s[ixs],evalf((2*x3)/(1+x3^2))):
        c3s := subs(x3=x3s[ixs],evalf((1-x3^2)/(1+x3^2))):
        s4s := subs(x4=x4s,evalf((2*x4)/(1+x4^2))):
        c4s := subs(x4=x4s,evalf((1-x4^2)/(1+x4^2))):
        s5s := subs(x5=x5s,evalf((2*x5)/(1+x5^2))):
        c5s := subs(x5=x5s,evalf((1-x5^2)/(1+x5^2))):
        S345:={s3=s3s,c3=c3s,s4=s4s,c4=c4s,s5=s5s,c5=c5s}:
        qs:=Transpose(Vt).DiagonalMatrix(map(s->1/s,S[1..8])).subs(S345,A8.p):

        s1s:=qs[5,1]:
        c1s:=qs[6,1]:
        s2s:=qs[7,1]:
        c2s:=qs[8,1]:
        S12:={s1=s1s,c1=c1s,s2=s2s,c2=c2s}:

        e61:=dhSimpl((<<1,0,0,0>|<0,-1,0,0>|<0,0,1,0>|<0,0,d2,1>>.
                        dhInv(iM22).M31.M32.M41.M42.M51.M52)[1..3,1..2],2):
        e62:=dhSimpl((<<1,0,0,0>|<0,-1,0,0>|<0,0,1,0>|<0,0,d2,1>>.
                        iM21.iM12.iM11.Mh.iM62.iM61)[1..3,1..2],2):
        e61v:=subs(S12,subs(S345,evalf(subs(MPh,e61)))):
        e62v:=subs(S12,subs(S345,evalf(subs(MPh,e62)))):
        e6:=convert(map(p->sort(p),e62v-e61v),Vector):
        t6:=[c6,s6,1]:
        M6:=PolyCoeffMatrix(M2L(e6),t6):
        cs6s:=MatrixInverse(-M6[1..2,1..2]).M6[1..2,3..3]:
        c6s:=cs6s[1,1]:
        s6s:=cs6s[2,1]:
        S6:={s6=cs6s[2,1],c6=cs6s[1,1]}:
        solutions := [op(solutions), evalf(subs(S12, S345, S6,
            <arctan(s1,c1), arctan(s2,c2), arctan(s3,c3),
             arctan(s4,c4), arctan(s5,c5), arctan(s6,c6)>))];
    end do;

    solutions:
end proc:
# Procedura pro prevod radianu na stupne
R2Deg := proc(x) evalf(x.180/Pi); end proc:
# Definice robotu, zadane polohy
# DH parametry robotu Mitsubishi Melfa RV6S
df := 0.01;
dfa := evalf(df.Pi/180);
alpha1:=-Pi/2+dfa: alpha2:=0+dfa: alpha3:=-Pi/2+dfa:
alpha4:=Pi/2+dfa: alpha5:=-Pi/2+dfa: alpha6:=0+dfa:
Mechanism := {
a1=85,  a2=280,  a3=100, a4=0+df, a5=0+df, a6=0+df,
d1=350, d2=0+df, d3=0+df, d4=315, d5=0+df, d6=85,
lambda1=cos(alpha1), mu1=sin(alpha1),
lambda2=cos(alpha2), mu2=sin(alpha2),
lambda3=cos(alpha3), mu3=sin(alpha3),
lambda4=cos(alpha4), mu4=sin(alpha4),
lambda5=cos(alpha5), mu5=sin(alpha5),
lambda6=cos(alpha6), mu6=sin(alpha6)};
# Parametry transformaci
# Pozice nulta a treti jsou stejne (t0 = t3, P0 = P3);
MhV0 := ImportMatrix("MhV0.mat", source=Matlab);
t0 := convert(ImportMatrix("t0.mat", source=Matlab), Vector):
MhV1 := ImportMatrix("MhV1.mat", source=Matlab);
t1 := convert(ImportMatrix("t1.mat", source=Matlab), Vector):
MhV2 := ImportMatrix("MhV2.mat", source=Matlab);
t2 := convert(ImportMatrix("t2.mat", source=Matlab), Vector):
# 
# Reseni ulohy
# Vypocet inverzni kinematicke ulohy
MhV := Matrix(4,4,[[lx,mx,nx,rx],[ly,my,ny,ry],[lz,mz,nz,rz],[0,0,0,1]]):
M0 := map(x->x[1]=x[2],convert(<convert(MhV[1..3,1..4],Vector)|convert(MhV0[1..3,1..4],Vector)>,listlist)):
M1 := map(x->x[1]=x[2],convert(<convert(MhV[1..3,1..4],Vector)|convert(MhV1[1..3,1..4],Vector)>,listlist)):
M2 := map(x->x[1]=x[2],convert(<convert(MhV[1..3,1..4],Vector)|convert(MhV2[1..3,1..4],Vector)>,listlist)):
solutions0 := IK6(Mechanism,M0);
solutions1 := IK6(Mechanism,M1);
solutions2 := IK6(Mechanism,M2);
# - s0 odpovida mym resenim inverzni kinematicke ulohy (vyjadrenym ve stupnich); predpokladam, ze jedno z nich by melo primo odpovidat hledanym offsetum
# - d1 a d2 odpovidaji rozdilum zadanych kloubovych souradnic a jim odpovidajicich, mnou vypoctenych, reseni inverzni kinematicke ulohy;
# - mezi nekterym ze sloupcu s0 a d1 (resp. s0 a d2) by melo dojit ke shode. tento sloupec odpovida hledanym offsetum
s0:=<R2Deg(solutions0[1])|R2Deg(solutions0[2])|R2Deg(solutions0[3])|R2Deg(solutions0[4])>;
d1:=<R2Deg(solutions1[1]) - t1|R2Deg(solutions1[2]) - t1|R2Deg(solutions1[3]) - t1|R2Deg(solutions1[4]) - t1>;
d2:=<R2Deg(solutions2[1]) - t2|R2Deg(solutions2[2]) - t2|R2Deg(solutions2[3]) - t2|R2Deg(solutions2[4]) - t2>;
# Procedura pro nalezeni sloupcu s nejvetsi mirou shody
findMinimalDifference := proc(s0, xx)
    minimum := infinity: minI := 0: minJ := 0:
    for i from 1 to Dimensions(s0)[2] do
        for j from 1 to Dimensions(d1)[2] do
            # norm of difference
            nod := Norm(s0[1..6,i] - xx[1..6,j]):
            if nod < minimum then
                minimum := nod:
                minI := i; minJ := j:
            end if:
        end do:
    end do:
    [minI, minJ, minimum]:
end proc:
# prvni dve cisla jsou indexy odpovidajici sloupcum, ktere se v s0 a d1 (resp. s0 a d2) nejvice shoduji; treti cislo je norma rozdilu techto sloupcu
d01 := findMinimalDifference(s0,d1);
d02 := findMinimalDifference(s0,d2);
d01[1];
# Nejlepsi nalezeny offset:
s0[1..6,d01[1]];
# spodni odhad "chyby" tohoto offsetu:
max(d01[3], d02[3]);
# 
# Zavery, pozorovani
# - Vyzkousel jsem vypocitat offsety pro hodnoty perturbaci z {1., 0.1, 0.01, 0.001, 0.0001}
# - Spodni odhad chyby offsetu byl nejmensi pro perturbaci df = 0.01
#?ExportMatrix
#ExportMatrix("sol0.mat", convert(solutions0, Matrix),target=Matlab);
#convert(solutions0, Matrix);

# 
