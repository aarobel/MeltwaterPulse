%% Params

clear all
close all
clc

global o
t_final = 1e5;

p.cross_section = 3.5e6; %the average length of the meridional ice sheet cross-section

p.h0 = 2.4;
p.h2 = 2.5;
p.Lc =1.4e6;

p.d0 = 0;
p.beta = 3e-3; %Depending on whether using PDD of snow or ice from Gregoir can be 0.005 to 0.015
p.s = 1.8e-3;    %Estimated by sight from Gregoire topography in supplement
% p.s2 = 1e-3;

p.da0dt = -0/1e5;
p.a00 = -1;    %Estimate by sight from Figure 3.13 in Ziemen Thesis 
p.ar = 1.4;    %Estimate by sight from Figure 3.13 in Ziemen Thesis

p.ic=[1600e3 1700e3]; %Initial conditions
p.tspan=[0,4e5]; %time ste

p.plot=0;
options = odeset('RelTol',1e-7,'AbsTol',1e-7,'NonNegative',1);
[time,Rout] = ode23tb(@(t,R) HMBF_2ODE_2sheets(t,R,p),p.tspan,p.ic,options);

R0 = Rout(end,:);
% figure;plot(time./1e3,Rout(:,1)./1e3,time./1e3,Rout(:,2)./1e3)

disp(['Vol is ' num2str(o.V(o.nt).*p.cross_section./1e9 /1e7) ' x 10^7 km^3'])
disp(['Height 1 is ' num2str(p.h0*sqrt(Rout(end,1))) ' m'])
disp(['Height 2 is ' num2str(p.h2*sqrt(Rout(end,2))) ' m'])


%% Steady-State
p.h0 = 2.4;
p.h2 = 2.6;
p.Lc =1.4e6;

p.d0 = 0;
p.beta = 3e-3; %Depending on whether using PDD of snow or ice from Gregoir can be 0.005 to 0.015
p.s = 1.8e-3;    %Estimated by sight from Gregoire topography in supplement
% p.s2 = 1e-3;

p.da0dt = -0/1e5;
p.a00 = -1;    %Estimate by sight from Figure 3.13 in Ziemen Thesis 
p.ar = 1.4;    %Estimate by sight from Figure 3.13 in Ziemen Thesis

%%Steady-state
p.da0dt = 0;
p.ic=R0;
p.tspan=linspace(0,2e5,10e3);

dt=p.tspan(2)-p.tspan(1);
Rout = p.ic;
for t = p.tspan

    rhs=HMBF_2ODE_2sheets(t,Rout,p);

    R1 = Rout(1) + dt*rhs(1);
    R2 = Rout(2) + dt*rhs(2);

    Rout(1) = R1;
    Rout(2) = R2;
end

R0 = Rout;


%% Vary forcing rate
clr = brewermap(5,'Dark2');q=0;
clr(3,:) = [0 0 0];
for da0dt = [-10 -7.5 -5 -2.5]./1e4
    q=q+1;
    p.da0dt = da0dt;
    p.tspan=linspace(0,2.5e4,2.5e3+1); %time steps
    p.ic = R0;
    Rout = R0;
    dt=p.tspan(2)-p.tspan(1);
    for t = p.tspan

        rhs=HMBF_2ODE_2sheets(t,Rout,p);

        R1 = Rout(1) + dt*rhs(1);
        R2 = Rout(2) + dt*rhs(2);

        Rout(1) = R1;
        Rout(2) = R2;
    end

    o.V(find(o.HS2(1:o.nt)==1,1))

    dVdt = diff(o.V(2:o.nt))./diff(o.t(2:o.nt));
    tm = (p.tspan(1:end-1)+p.tspan(2:end))./2;
    a0 = o.a0(1:o.nt);
    a0m = (a0(1:end-1)+a0(2:end))./2;
    
    figure(5);subplot(2,1,1)
    plot(tm(2:end)./1e3,-dVdt(2:end).*p.cross_section./1e9./4.6e3,'-','linewidth',5,'Color',clr(q,:));hold on
    
    subplot(2,1,2)
    plot(-tm(2:end).*p.da0dt,dVdt(2:end).*p.cross_section./1e9./4.6e3./a0m(3:end),'-','linewidth',5,'Color',clr(q,:));hold on
end

figure(5);set(5,'units','pixels','position',[0 0 1002 802]);
subplot(2,1,1)
ylabel('Ice loss rate (SLE; cm/yr)','fontsize',18)
set(gca,'fontsize',24);xlabel('Time (kyr)','fontsize',18);
xlim([0 12]);ylim([0 3])
text(0.01,0.99,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)
drawnow

subplot(2,1,2)
ylabel('Loss rate norm. by forcing (dV/dt / a_0(t))','fontsize',18)
set(gca,'fontsize',24);xlabel('Time scaled by forcing rate (t \times a_{tr})','fontsize',18);
xlim([0 7]);%ylim([0 3.5])
text(0.01,0.99,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)

%% Quasi-equilibrium run
q=q+1;
p.da0dt = -0.01/1e4;
p.tspan=linspace(0,500e4,500e3+1); %time steps
p.ic = R0;
Rout = R0;
dt=p.tspan(2)-p.tspan(1);
for t = p.tspan

    rhs=HMBF_2ODE_2sheets(t,Rout,p);

    R1 = Rout(1) + dt*rhs(1);
    R2 = Rout(2) + dt*rhs(2);

    Rout(1) = R1;
    Rout(2) = R2;
end


dVdt = diff(o.V(2:o.nt))./diff(o.t(2:o.nt));
tm = (p.tspan(1:end-1)+p.tspan(2:end))./2;
a0 = o.a0(1:o.nt);
a0m = (a0(1:end-1)+a0(2:end))./2;
    
figure(5)

h1=subplot(2,1,1);
plot(tm(2:end)./1e3,-dVdt(2:end).*p.cross_section./1e9./4.6e3,'b-','linewidth',5);hold on
xlim([0 1811]);breakxaxis([12 1805])
set(h1,'xtick',[0:2:10,1810])

subplot(2,1,2)
plot(-tm(2:end).*p.da0dt,dVdt(2:end).*p.cross_section./1e9./4.6e3./a0m(3:end),'b-','linewidth',5);hold on
legend('a_{tr} = 10 \times 10^{-4} m/yr^2','a_{tr} = 7.5 \times 10^{-4} m/yr^2','a_{tr} = 5 \times 10^{-4} m/yr^2','a_{tr} = 2.5 \times 10^{-4} m/yr^2','a_{tr} = 1 \times 10^{-6} m/yr^2','Location','NorthEast')
