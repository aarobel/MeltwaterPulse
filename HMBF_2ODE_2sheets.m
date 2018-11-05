function rhs=HMBF_2ODE_2sheets(t,R,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% right hand side of evolution equations

persistent nt
global o

%% Climate warming
p.a0 = p.a00 + p.da0dt*t;

%% Read in prognostic variables
R1=max([R(1),0]);
R2=max([R(2),0]);

%% Calulate intersection point, runoff lines and corresponding heavyside functions

alpha = 1/(1 + (p.h0/p.h2)^2);
x_i = alpha*((p.h0/p.h2)^2 * R1 - R2 + p.Lc);
% x_i = max([x_i 0]);

HS = heaviside(R1+R2-p.Lc); %1 when sheets intersect, 0 when separated
if(R1 < 0.01*p.Lc || R2 < 0.01*p.Lc)
    HS = 0;
end

x_i12 = x_i*HS + R1*(1-HS); %intersection points (saddle)
x_i34 = (p.Lc-x_i)*HS + R2*(1-HS); 

x_r1 = (p.h0/(p.s^2))*sqrt(((p.ar-p.a0)/p.beta - p.d0)*p.s + ((p.h0^2) / 4) + R1*(p.s^2)) -...
    ((p.ar-p.a0)/p.beta - p.d0)/p.s - 0.5*(p.h0/p.s)^2; %location of saturation

x_r2_p = R1 - (1/p.h0^2)*((p.ar-p.a0)/p.beta - p.d0)^2;

HS2 = heaviside(x_i12-x_r2_p);  %0 when ice sheets have no interior melt, 1 when there is
x_r2 = x_r2_p*HS2 + x_i12*(1-HS2);

x_r3_p = R2 - (1/p.h2^2)*((p.ar-p.a0)/p.beta - p.d0)^2; 
x_r3 = x_r3_p*HS2 + x_i34*(1-HS2);

x_r4 = (p.h2/(p.s^2))*sqrt(((p.ar-p.a0)/p.beta - p.d0)*p.s + ((p.h2^2)/4) + R2*(p.s^2)) -...
    ((p.ar-p.a0)/p.beta - p.d0)/p.s - 0.5*(p.h2/p.s)^2;

%ensure that location of saturation maxes out when runoff line is above ice
%sheet max (thanks to reviewer #2)
x_r1 = max([0 x_r1]);
x_r2 = max([0 x_r2]);
x_r3 = max([0 x_r3]);
x_r4 = max([0 x_r4]);

%% Calculate SMB on both ice sheets
B1 = p.ar*x_r1 + p.a0*(R1-x_r1) + p.d0*p.beta*(R1-x_r1) - 0.5*p.beta*p.s*(R1^2 - x_r1^2) + ...
    (2/3)*p.h0*p.beta*(R1-x_r1)^(3/2);
B2 = (p.ar-p.a0)*x_r2 + p.a0*x_i12 + (2/3)*p.beta*p.h0*((R1-x_r2)^(3/2) -...
    (R1-x_i12)^(3/2)) + p.beta*p.d0*(x_i12-x_r2);

B3 = (p.ar-p.a0)*x_r3 + p.a0*x_i34 + (2/3)*p.beta*p.h2*((R2-x_r3)^(3/2) -...
    (R2-x_i34)^(3/2)) + p.beta*p.d0*(x_i34-x_r3);
B4 = x_r4*(p.ar-p.a0) + p.a0*R2 + p.d0*p.beta*(R2 - x_r4) - 0.5*p.beta*p.s*(R2^2 - x_r4^2) +...
    (2/3)*p.h2*p.beta*(R2 - x_r4)^(3/2);

B1B2 = B1+B2;
B3B4 = B3+B4;   

%% Calculate volume integrals
V1 = (2/3)*p.h0*(R1^(3/2));
V2 = (2/3)*p.h0*(R1^(3/2)) - (2/3)*p.h0*((R1-x_i12)^(3/2));
V3 = (2/3)*p.h2*(R2^(3/2));
V4 = (2/3)*p.h2*(R2^(3/2)) - (2/3)*p.h2*((R2-x_i34)^(3/2));

%% Calculate coefficients for terms
a = 2*p.h0*sqrt(R1) - p.h0*(alpha^(3/2))*sqrt(R1+R2-p.Lc);
b = p.h0*(alpha^(3/2))*sqrt(R1+R2-p.Lc);
c = 2*p.h2*sqrt(R2) - p.h2*((1-alpha)^(3/2))*sqrt(R1+R2-p.Lc);
d = p.h2*((1-alpha)^(3/2))*sqrt(R1+R2-p.Lc);

%% Calculate evolution equations
if(HS==1)
    dRdt(1) = (c*B1B2 + b*B3B4)/(a*c - b*d);
    dRdt(2) = (d*B1B2 + a*B3B4)/(a*c - b*d);
else
    dRdt(1) = B1B2/(2*p.h0*sqrt(R1));
    dRdt(2) = B3B4/(2*p.h2*sqrt(R2));
end

if(R1<=0)
    dRdt(1) = 0;
end

if(R2<=0)
    dRdt(2) = 0;
end

%% save diagnostics to output variable:
% first, advance time step used for saving diagnostics.
if isempty(nt); nt=1; else; nt=nt+1; end
if nt==1; init=1; else; init=0; end

max_nt=10000;

%% matlab sometimes makes a trial step, and then goes back to get
%% more accurate results, take care of this by not saving trial
%% steps in output arrays:
%%fprintf(1,'%g, %d; ',t, nt);
while nt>2 && t<=o.t(nt-1);
  nt=nt-1;
  %%fprintf(1,'\n !! \n');
end
o.nt=nt;


%% Save stuff
if init; ZZ=zeros(max_nt,1)*NaN; end
if init; o.R1=ZZ; end;   o.R1(nt) = R1;
if init; o.R2=ZZ; end;   o.R2(nt) = R2;

if init; o.B1=ZZ; end;   o.B1(nt) = B1;
if init; o.B2=ZZ; end;   o.B2(nt) = B2;
if init; o.B3=ZZ; end;   o.B3(nt) = B3;
if init; o.B4=ZZ; end;   o.B4(nt) = B4;

if init; o.HS1=ZZ; end;   o.HS1(nt) = HS;
if init; o.HS2=ZZ; end;   o.HS2(nt) = HS2;
if init; o.xi1=ZZ; end;   o.xi1(nt) = x_i12;
% if init; o.HS3=ZZ; end;   o.HS3(nt) = HS3;

if init; o.a0=ZZ; end;   o.a0(nt) = p.a0;

if init; o.V1=ZZ; end;   o.V1(nt) = V1;
if init; o.V2=ZZ; end;   o.V2(nt) = V2;
if init; o.V3=ZZ; end;   o.V3(nt) = V3;
if init; o.V4=ZZ; end;   o.V4(nt) = V4;

if init; o.V=ZZ; end;   o.V(nt) = V1+V2+V3+V4;
if init; o.t=ZZ; end;   o.t(nt)=t;
rhs= dRdt';
