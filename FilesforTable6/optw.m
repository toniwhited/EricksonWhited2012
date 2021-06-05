%subroutine to define the f vector for the optimal weighting matrix.


function[f] = optw(a,mmm,ooo,zzz,n,neq,estim,wf);
f=zeros(n,neq);


f(:,1)=mmm{1}-ooo{1};
f(:,2)=mmm{2}-ooo{2};
f(:,3)=mmm{3}-ooo{3};
f(:,4)=mmm{4}-ooo{4};
f(:,5)=mmm{5}-ooo{5};
if estim>1;
f(:,6)=mmm{6}-ooo{6};
f(:,7)=mmm{7}-ooo{7};
f(:,8)=mmm{8}-ooo{8};
end;
if estim>2;
f(:,9)=mmm{9}-ooo{9};
f(:,10)=mmm{10}-ooo{10};
f(:,11)=mmm{11}-ooo{11};
f(:,12)=mmm{12}-ooo{12};
f(:,13)=mmm{13}-ooo{13};
f(:,14)=mmm{14}-ooo{14};
end;
if estim>3;
f(:,15)=mmm{15}-ooo{15};
f(:,16)=mmm{16}-ooo{16};
f(:,17)=mmm{17}-ooo{17};
f(:,18)=mmm{18}-ooo{18};
f(:,19)=mmm{19}-ooo{19};
f(:,20)=mmm{20}-ooo{20};
f(:,21)=mmm{21}-ooo{21};
end;
if estim>4;
f(:,22)=mmm{22}-ooo{22};
f(:,23)=mmm{23}-ooo{23};
f(:,24)=mmm{24}-ooo{24};
f(:,25)=mmm{25}-ooo{25};
f(:,26)=mmm{26}-ooo{26};
f(:,27)=mmm{27}-ooo{27};
f(:,28)=mmm{28}-ooo{28};
f(:,29)=mmm{29}-ooo{29};
end;

%this part makes the standard error adjustment.

if wf==2;


inezz    = zzz{1};
zy_d     = zzz{2};
zx_d     = zzz{3};
eyx_z    = zzz{4};
ey2_z    = zzz{5};
ex2_z    = zzz{6};
ey3_z    = zzz{7};
ey2x_z   = zzz{8};
eyx2_z   = zzz{9};
ex3_z    = zzz{10};
ey4_z    = zzz{11};
ey3x_z   = zzz{12};
ey2x2_z  = zzz{13};
eyx3_z   = zzz{14};
ex4_z    = zzz{15};
ey5_z    = zzz{16};
ey4x_z   = zzz{17};
ey3x2_z  = zzz{18};
ey2x3_z  = zzz{19};
eyx4_z   = zzz{20};
ex5_z    = zzz{21};
ey6_z    = zzz{22};
ey5x_z   = zzz{23};
ey4x2_z  = zzz{24};
ey3x3_z  = zzz{25};
ey2x4_z  = zzz{26};
eyx5_z   = zzz{27};
ex6_z    = zzz{28};

ei=zeros(n,neq);

ei(:,4)=(-2*eyx_z'*inezz*(zy_d')-ey2_z'*inezz*(zx_d'))';
ei(:,5)=(-ex2_z'*inezz*zy_d'-2*eyx_z'*inezz*zx_d')';


if estim>1;
ei(:,6)=(-3*ey2x_z'*inezz*zy_d'-ey3_z'*inezz*zx_d')';
ei(:,7)=(-2*eyx2_z'*inezz*zy_d'-2*ey2x_z'*inezz*zx_d')';
ei(:,8)=(-ex3_z'*inezz*zy_d'-3*eyx2_z'*inezz*zx_d')';
end;


if estim>2;
ei(:,9)=(-3*ex2_z'*inezz*zx_d')';
ei(:,10)=(-3*ey2_z'*inezz*zy_d')';
ei(:,11)=(-4*ey3x_z'*inezz*zy_d'-ey4_z'*inezz*zx_d')';
ei(:,12)=(-3*ey2x2_z'*inezz*zy_d'-2*ey3x_z'*inezz*zx_d')';
ei(:,13)=(-2*eyx3_z'*inezz*zy_d'-3*ey2x2_z'*inezz*zx_d')';
ei(:,14)=(-ex4_z'*inezz*zy_d'-4*eyx3_z'*inezz*zx_d')';
end;


if estim>3;
ei(:,15)=(-4*ey3_z'*inezz*zy_d')';
ei(:,16)=(-4*ex3_z'*inezz*zx_d')';


ei(:,17)=(-5*ey4x_z'*inezz*zy_d'-ey5_z'*inezz*zx_d')';
ei(:,18)=(-4*ey3x2_z'*inezz*zy_d'-2*ey4x_z'*inezz*zx_d')';
ei(:,19)=(-3*ey2x3_z'*inezz*zy_d'-3*ey3x2_z'*inezz*zx_d')';
ei(:,20)=(-2*eyx4_z'*inezz*zy_d'-4*ey2x3_z'*inezz*zx_d')';
ei(:,21)=(-ex5_z'*inezz*zy_d'-5*eyx4_z'*inezz*zx_d')';
end;


if estim>4;
ei(:,22)=(-5*ey4_z'*inezz*zy_d')';
ei(:,23)=(-5*ex4_z'*inezz*zx_d')';
ei(:,24)=(-6*ey5x_z'*inezz*zy_d'-ey6_z'*inezz*zx_d')';
ei(:,25)=(-5*ey4x2_z'*inezz*zy_d'-2*ey5x_z'*inezz*zx_d')';
ei(:,26)=(-4*ey3x3_z'*inezz*zy_d'-3*ey4x2_z'*inezz*zx_d')';
ei(:,27)=(-3*ey2x4_z'*inezz*zy_d'-4*ey3x3_z'*inezz*zx_d')';
ei(:,28)=(-2*eyx5_z'*inezz*zy_d'-5*ey2x4_z'*inezz*zx_d')';
ei(:,29)=(-ex6_z'*inezz*zy_d'-6*eyx5_z'*inezz*zx_d')';
end;


f=f+ei;
end;
