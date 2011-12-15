# IRO-2011-IK-3-axes-6-dof-solution-assignment.mws - Inverse Kinematics of  a manipulator with three axes of motion
# T. Pajdla, 4 Dec 2011
# Assignment
#  Packages & settings
restart;
with(LinearAlgebra):
with(PolynomialTools):
interface(rtablesize=50):
Digits:=30:
#  DH-Kinematics
# Functions for DH kinematics
# Joint transformations:
# Two one-parametric motions transformatin in DH-convention(phi, theta, a,d) indexed by i
# c = cos(phi), s = sin(phi), P = cos(alpha), R = sin(alpha)
dhTs := proc(i)
local M1, M2;
   M1:=Matrix(4,4,[[+cat(`c`,i),-cat(`s`,i),0,         0],
                  [ +cat(`s`,i),+cat(`c`,i),0,         0],
                  [           0,          0,1,cat(`D`,i)],
                  [           0,          0,0,         1]]);
   M2:=Matrix(4,4,[[1,        0,          0,cat(`A`,i)],  
                  [ 0,+cat(`P`,i),-cat(`R`,i),         0],
                  [ 0,+cat(`R`,i),+cat(`P`,i),         0],
                  [ 0,          0,          0,         1]]);
   [M1,M2];
end proc:
# Inverse of the DH-convention for one-aprametric DH rigid motion transformations
dhInvs := proc(M)
   local M1, M2;
   M1 := M[1];
   M2 := M[2];
   [simplify(MatrixInverse(M2),{M2[3,2]^2+M2[3,3]^2=1}),
    simplify(MatrixInverse(M1),{M1[1,1]^2+M1[2,1]^2=1})];
end proc:
# Rigid motion transformatin in DH-convention(phi, theta, a,d) indexed by i
# c = cos(phi), s = sin(phi), P = cos(alpha), R = sin(alpha)
dhT := proc(i)
local M;
M:=Matrix(4,4,[[+cat(`c`,i),-cat(`s`,i)*cat(`P`,i),+cat(`s`,i)*cat(`R`,i),+cat(`A`,i)*cat(`c`,i)],
               [+cat(`s`,i),+cat(`c`,i)*cat(`P`,i),-cat(`c`,i)*cat(`R`,i),+cat(`A`,i)*cat(`s`,i)],
               [         0,             cat(`R`,i),            cat(`P`,i),            cat(`D`,i)],
               [         0,                      0,                     0,                     1]]);
end proc:
# Inverse of the DH-convention rigid motion transformation
dhInv := proc(M)
   simplify(MatrixInverse(M),{M[1,1]^2+M[2,1]^2=1,M[3,2]^2+M[3,3]^2=1});
end proc:
# Simplify using the trigonometric indentities c^2+s^2=1 & P^2+R^2=1
dhSimpl := proc(M,i)
   simplify(M,{cat(`c`,i)^2+cat(`s`,i)^2=1,cat(`P`,i)^2+cat(`R`,i)^2=1});
end proc:
# Simplify a general motion matrix using rotation matrix identities in columns
rcSimpl := proc(M,R)
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
# Simplify a general motion matrix using rotation matrix identities in rows
rrSimpl := proc(M,R)
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
# Matrix representation of a set of polynomials 
PolyCoeffMatrix:=proc(S,m) 
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
# n x 1 matrix to a list conversion
M2L:=proc(M) 
        convert(convert(M,Vector),list);
end proc:
# Highlit non-zero entries
spy:=proc(A)
   map(x->`if`(simplify(x)=0,`.`, `if`(simplify(x)=1,1,`*`)) ,A):
end proc:
#
## Monomials of a set of polynomials
#
PolyVarMonomials:=proc(S::list,Ord) 
local v,m,i,c,q;
        v:={op(Ord)};
        m:=[];
        for i from 1 to nops(S) do
                c:=[coeffs(expand(S[i]),v,'q')];
                m:=[op(m),q];
        end do;
        m:=ListTools[MakeUnique](m);
        sort(m,(t1,t2)->Groebner[TestOrder](t2,t1,Ord));
end proc:
#
#  IK solution for a general case
# The formulation:
# 
# (1) M1*M2*M3 = Mh
# 
# or equivalently
# 
# (2) M1*M2 = Mh*M3^{-1}

# i.e.
# 
# (3) M11*M12*M21*M22 = Mh*M32^{-1}*M31^{-1}
# 
# Construct motion matrices:
M11:=dhTs(1)[1]:
M12:=dhTs(1)[2]:
M21:=dhTs(2)[1]:
M22:=dhTs(2)[2]:
iM31:=dhInvs(dhTs(3))[2]:
iM32:=dhInvs(dhTs(3))[1]:
Mh:=Matrix(4,4,[[Lx,Mx,Nx,Qx],[Ly,My,Ny,Qy],[Lz,Mz,Nz,Qz],[0,0,0,1]]):
Lhs:=M11.M12.M21.M22:
Rhs:=Mh.iM32.iM31:
# See the equation (3) in separate matrices
M11,M12,M21,M22, "=", Mh,iM32,iM31;
# as well as in their  product:
Lhs, "=", Rhs;
# Let's get all 12 equations.
lh:=convert(Lhs[1..3,1..4],Vector):
rh:=convert(Rhs[1..3,1..4],Vector):
lh, "=", rh;

# Add trigonometric identities
lh1:=<lh, c1^2+s1^2, c2^2+s2^2, c3^2+s3^2>:
rh1:=<rh, 1, 1, 1>:
# Move all to the left-hand side.
E6:= lh1 - rh1:
E6, "=", ZeroMatrix(12,1);
# The above equations describe our mechanism completely.
# 
# Computing the general solution of this problem can't be anymore done in a symbolic form
# bacause equations take quite different form for general and special values of parameters
# of the mechanism. Therefore, we will design a symbolic-numerical method.
# 
#  Symbolic-numerical solution
# Simulate a mechanism.
# (4)         M1*M2*M3 = Mh
M1x := dhT(1):
M2x := dhT(2):
M3x := dhT(3):
Mhx := M1x.M2x.M3x:
"Mhx =", M1x, M2x, M3x;
     
# Set the fixed parameters of the mechanism
alphasMech := {al1=-Pi/3,al2=Pi/5,al3=-Pi/7}; 
Mech := {A1=1,A2=1,A3=1,P1=cos(al1),R1=sin(al1),P2=cos(al2),R2=sin(al2),P3=cos(al3),R3=sin(al3)};  
Mech:=map(f->subs(alphasMech,f),Mech);
# Set the controled parameters of the mechanism
thetasPos := {th1=0.1,th2=0.2,th3=0.3}; # Ad 0), 1)
#thetasPos := {th1=0.0,th2=0.2,th3=0.0}; # Ad 2)
#thetasPos := {th1=0.1,th2=0.0,th3=0.3}; # Ad 3)
Pos := {D1=1,D2=2,D3=3,c1=cos(th1),s1=sin(th1),c2=cos(th2),s2=sin(th2),c3=cos(th3),s3=sin(th3)};  
Pos:=map(f->subs(thetasPos,f),Pos);

# Solve  DK  of this mechanism.
MhX:=evalf(subs(Pos,subs(Mech,Mhx)));

# Extract the values from MhX in a suitable form for further substitutions.
MhS:={Lx=MhX[1,1],Mx=MhX[1,2],Nx=MhX[1,3],Qx=MhX[1,4],
      Ly=MhX[2,1],My=MhX[2,2],Ny=MhX[2,3],Qy=MhX[2,4],
      Lz=MhX[3,1],Mz=MhX[3,2],Nz=MhX[3,3],Qz=MhX[3,4]}:

# Copy the general equations derived above to a new matrix E
E := Matrix(15, 1,[[c1*c2-s1*P1*s2-Lx*c3+(Mx*P3-Nx*R3)*s3],
                      [s1*c2+c1*P1*s2-Ly*c3+(My*P3-Ny*R3)*s3],
                      [R1*s2-Lz*c3+(Mz*P3-Nz*R3)*s3],
                      [(-c1*s2-s1*P1*c2)*P2+s1*R1*R2-Lx*s3-(Mx*P3-Nx*R3)*c3],
                      [(-s1*s2+c1*P1*c2)*P2-c1*R1*R2-Ly*s3-(My*P3-Ny*R3)*c3],
                      [R1*c2*P2+P1*R2-Lz*s3-(Mz*P3-Nz*R3)*c3],
                      [-(-c1*s2-s1*P1*c2)*R2+s1*R1*P2-Mx*R3-Nx*P3],
                      [-(-s1*s2+c1*P1*c2)*R2-c1*R1*P2-My*R3-Ny*P3],
                      [-R1*c2*R2+P1*P2-Mz*R3-Nz*P3],
                      [(c1*c2-s1*P1*s2)*A2+s1*R1*D2+c1*A1+(Mx*R3+Nx*P3)*D3+Lx*A3-Qx],
                      [(s1*c2+c1*P1*s2)*A2-c1*R1*D2+s1*A1+(My*R3+Ny*P3)*D3+Ly*A3-Qy],
                      [R1*s2*A2+P1*D2+D1+(Mz*R3+Nz*P3)*D3+Lz*A3-Qz],
                      [c1^2+s1^2-1],
                      [c2^2+s2^2-1],
                      [c3^2+s3^2-1]]);
# Substitute mechanism and position values into the equations 
E1:=simplify(subs(MhS,subs(Mech,E)));
# Solve for c2
# Extract monomials in such an order that c2 is the last variable
t1:=<<PolyVarMonomials(M2L(E1),plex(D3,D2,D1,s3,c3,s1,c1,s2,c2))>>:
Transpose(t1):

# Express the equations in the matrix form w.r.t. the monomials
T1:=PolyCoeffMatrix(M2L(E1),M2L(t1)):
<Transpose(t1),spy(T1)>:

# Simplify the equations by using Linear Algebra
T1:=GaussianElimination(T1):
<Transpose(t1),spy(T1)>;

Dimension(T1);

# tady Rank(T1) je jenom pomucka pro nas, melo by tady byt cislo radky...
# The last column is a polynomial in 1 variable, i.e. c2, solve it
Rank(T1);
T1[15,21];
eq1:=convert(T1[Rank(T1),1..ColumnDimension(T1)].t1,set);
#eq1:=convert(T1[15,1..ColumnDimension(T1)].t1,set);
c2s:=solve(eq1);
# Solve for s2 
Dimension(E1);
E1[1];
E2[1];
# Substitute the solution for s2 back to equations and repeat the above procedure
E2:=simplify(subs(c2s[1],T1.t1)):
t2:=<<PolyVarMonomials(M2L(E2),plex(D3,D2,D1,s3,c3,s1,c1,s2))>>:
T2:=PolyCoeffMatrix(M2L(E2),M2L(t2)):
# <Transpose(t2),spy(T2)>;
# T2[14,18];
T2:=GaussianElimination(T2):
<Transpose(t2),spy(T2)>;
# Dimension(T2);
# Rank(T2);
T2[13,18]; # v T2 na predposlednim a predpredposlednim radku jsou to skoro nuly
# The last column is a polynomial in 1 variable, i.e. s2, solve it
eq2:=convert(T2[Rank(T2),1..ColumnDimension(T2)].t2,set);
#eq2:=convert(T2[13,1..ColumnDimension(T2)].t2,set);
s2s:=solve(eq2);
# Solve for c1
E3:=simplify(subs(s2s[1],T2.t2)):   # nezapomenout taky na druhe reseni... pouzit i s2s[2]
t3:=<<PolyVarMonomials(M2L(E3),plex(D3,D2,D1,s3,c3,s1,c1))>>:
T3:=PolyCoeffMatrix(M2L(E3),M2L(t3)):
    #<Transpose(t3),spy(T3)>;
T3:=GaussianElimination(T3):
    <Transpose(t3),spy(T3)>;
`Dim`,`=`, Dimension(T3);
`R`,`=`, Rank(T3);
    T3[9,13];
eq3:=convert(T3[Rank(T3),1..ColumnDimension(T3)].t3,set);
#eq3:=convert(T3[9,1..ColumnDimension(T3)].t3,set);
c1s:=solve(eq3);
# Solve for s1
E4:=simplify(subs(c1s[1],T3.t3)):
t4:=<<PolyVarMonomials(M2L(E4),plex(D3,D2,D1,s3,c3,s1))>>:
T4:=PolyCoeffMatrix(M2L(E4),M2L(t4)):
    #<Transpose(t4),spy(T4)>;
T4:=GaussianElimination(T4):
    <Transpose(t4),spy(T4)>;
`Dim`,`=`, Dimension(T4);
`R`,`=`, Rank(T4);
    T4[8,11];
eq4:=convert(T4[Rank(T4),1..ColumnDimension(T4)].t4,set);
#eq4:=convert(T4[8,1..ColumnDimension(T4)].t4,set);
s1s:=solve(eq4);
# Solve for c3
E5:=simplify(subs(s1s[1],T4.t4)):
t5:=<<PolyVarMonomials(M2L(E5),plex(D3,D2,D1,s3,c3))>>:
T5:=PolyCoeffMatrix(M2L(E5),M2L(t5)):
    #<Transpose(t5),spy(T5)>;
T5:=GaussianElimination(T5):
    <Transpose(t5),spy(T5)>;
`Dim`,`=`, Dimension(T5);
`R`,`=`, Rank(T5);
    #T5[9,11];
eq5:=convert(T5[Rank(T5),1..ColumnDimension(T5)].t5,set);
#eq5:=convert(T5[6,1..ColumnDimension(T5)].t5,set);
c3s:=solve(eq5);
# Solve for s3
E6:=simplify(subs(c3s[1],T5.t5)):
t6:=<<PolyVarMonomials(M2L(E6),plex(D3,D2,D1,s3))>>:
T6:=PolyCoeffMatrix(M2L(E6),M2L(t6)):
    #<Transpose(t6),spy(T6)>;
T6:=GaussianElimination(T6):
    <Transpose(t6),spy(T6)>;
`Dim`,`=`, Dimension(T6);
`R`,`=`, Rank(T6);
    T6[5,6];
eq6:=convert(T6[Rank(T6),1..ColumnDimension(T6)].t6,set);
#eq6:=convert(T6[5,1..ColumnDimension(T6)].t6,set);
s3s:=solve(eq6);
# Solve for D1, D2, D3
E7:=simplify(subs(s3s[1],T6.t6)):
t7:=<<PolyVarMonomials(M2L(E7),plex(D3,D2,D1))>>:
T7:=PolyCoeffMatrix(M2L(E7),M2L(t7)):
    #<Transpose(t7),spy(T7)>;
Te7:=ReducedRowEchelonForm(T7[1..3,1..4]):
    <Transpose(t7),spy(T7)>;
    #Te7[1..3,1..4];
D1s:=-Te7[3,4];
D2s:=-Te7[2,4];
D3s:=-Te7[1,4];
# Compare the results
evalf(subs(s3s, c3s, s2s, c2s, s1s, c1s, [D1s, D2s, D3s, arctan(s3,c3), arctan(s2,c2), arctan(s1,c1)]));
s2s;
errors := subs(evalf(Pos),[D3,D2,D1,s3,c3,s2,c2,s1,c1])-subs({D3=D3s,D2=D2s,D1=D1s} union s3s union c3s union s2s[1] union c2s union s1s union c1s,[D3,D2,D1,s3,c3,s2,c2,s1,c1]);
# Uloha mela dve ruzna reseni. Puvodnim hodnotam odpovida prvni reseni, odvozene z prvniho reseni pro s2 (s2s[1]).
# So it works!
# Now try:
# 
# 1)  Digits:=20, 10, 5,2
# 2)  {th1=0.0,th2=0.2,th3=0.0}
# 3)  {th1=0.1,th2=0.0,th3=0.3}
# 
# Was there any problem? If so, why and what can be done?
# Odpovedi
# 2. a. podle zadani na webu 
# al1=-Pi/3	al2=Pi/5	al3=-Pi/7	A1=1	A2=1	A3=1
# th1=0.1	th2=0.2	th3=0.3	D1=1	D2=2	D3=3
# vliv volby Digits na presnost ziskanych reseni
# Vypoctena byla dve reseni. Zadane poloze odpovida prvni z nich. Odchylky tohoto prvniho reseni od zadaneho jsou pro ruzne hodnoty promenne Digits uvedny v sekci 3.
# 2. b
# 
# al1=-Pi/3	al2=Pi/5	al3=-Pi/7	A1=1	A2=1	A3=1
# th1=0.1	th2=0.0	th3=0.3	D1=1	D2=2	D3=3
# 
# Typesetting:-Parse:-ConvertTo1D, "first argument to _Inert_ASSIGN must be assignable";
# 3. Vliv Digits na presnost ziskanych vysledku
# Ziskane hodnoty errors pro zadani 2.a.:
# - pro Digits = 30:     errors := [3.8*10^(-28), -1.9*10^(-28), -2.3*10^(-28), -1.3*10^(-29), 0., 5.*10^(-30), -1.*10^(-30), 7.8*10^(-30), -2.*10^(-30)];
# - pro Digits = 20:     errors := [5.*10^(-19), -7.*10^(-19), -2.*10^(-19), -1.1*10^(-19), 4.*10^(-20), 3.*10^(-20), -1.*10^(-20), 7.3*10^(-20), 4.*10^(-20)];
# - pro Digits = 10:     errors := [3.6*10^(-8), -2.2*10^(-8), -2.0*10^(-8), -3.*10^(-10), 1.5*10^(-9), 3.*10^(-10), -1.*10^(-10), 9.4*10^(-10), 0.];
# - pro Digits = 5:       errors := [-0.14e-2, 0.9e-3, 0.83e-3, 0.10e-3, 0.3e-4, -0.3e-4, 0.1e-4, -0.60e-4, 0.4e-4];
# - pro Digits = 2:       errors := [-6.6, 1.0, 6.3, -.51, 1.1, .88, 0., .70, .51];

# Je zrejme, ze zvolena hodnota promenne Digits ovlivnuje presnost ziskanych reseni.
# 
# Typesetting:-Parse:-ConvertTo1D, "first argument to _Inert_ASSIGN must be assignable";
# 
# {th1=0.0,th2=0.2,th3=0.0}
# Vypada to, ze s touto alternativou neni problem, jedine jsem si vsimnul,, ze pro s3 a s1 vysel rozdil vysledku oproti puvodni hodnote na vyssi pocet platnych mist:
#     errors := [-1.1*10^(-28), 8.*10^(-29), 6.2*10^(-29), 1.12237272943531195877167605453*10^(-29), 0., -5.*10^(-30), 1.*10^(-30), -8.28206014728042101468891029684*10^(-30), 0.];
# {th1=0.1,th2=0.0,th3=0.3}
# 

