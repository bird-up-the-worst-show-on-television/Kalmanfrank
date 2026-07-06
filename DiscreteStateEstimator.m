clear all
syms Vy Vx r Psi Ax Rbias Abias Vypri Vxpri rpri Psipri Axpri Rbiaspri Abiaspri Mv B D C B E Iz delta a b Ts 

AlphaF = delta - atan((-Vypri+a*(rpri+Rbiaspri))/Vxpri);
AlphaR =       -atan((-Vypri-b*(rpri+Rbiaspri))/Vxpri);


Ff = D*sin(C*atan(B*AlphaF-E*(B*AlphaF-atan(B*AlphaF))));
Fr = D*sin(C*atan(B*AlphaR-E*(B*AlphaR-atan(B*AlphaR))));
 
Vy    = Vypri  + ((Fr+Ff*cos(delta))/Mv-Vxpri*rpri)*Ts;
Vx    = Vxpri  + ((Axpri+Abiaspri)-Vypri*rpri)*Ts; 
r     = (rpri)   + ((Ff*cos(delta)*a-Fr*b)*Ts)/Iz;
Psi   = Psipri + (rpri)*Ts;
Ax    = Axpri+Abiaspri;
Rbias = Rbiaspri;
Abias = Abiaspri;

Xk = [Vy;Vx;r;Psi;Ax;Rbias;Abias];
Fk = jacobian(Xk,[Vypri,Vxpri,rpri,Psipri,Axpri,Rbiaspri,Abiaspri])

matlabFunction(Xk,   'File','scripts/LittleF');
matlabFunction(Fk,   'File','scripts/BigF');

