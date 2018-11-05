%% Steady State

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


%% h0 vs h2 2D
np = 50;
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
q=0;
n=0;

for h0 = linspace(2.1,2.6,np)
    p.h0 = h0;
    q=q+1; 
    n=0;
    q
    for h2 = linspace(2.3,3,np)
        p.h2 = h2;
        n=n+1;
        
        h0s(q,n) = h0;
        h2s(q,n) = h2;
        
        
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
        
        if(R1<=0 || R2<=0)
            if(R1<R2)
                R1ss(q,n) = nan;
                R2ss(q,n) = Rout(end,2);
                Vss(q,n) = nan;%(o.V3(o.nt)+o.V4(o.nt)).*p.cross_section./1e9;
                DVmax(q,n) = nan;
                DVtime(q,n) = nan;
                MWPdur(q,n) = nan;
                continue
            else
                R1ss(q,n) = Rout(end,1);
                R2ss(q,n) = nan;
                Vss(q,n) = nan;%(o.V1(o.nt)+o.V2(o.nt)).*p.cross_section./1e9;
                DVmax(q,n) = nan;
                DVtime(q,n) = nan;
                MWPdur(q,n) = nan;
                continue
            end
        else
            R1ss(q,n) = Rout(end,1);
            R2ss(q,n) = Rout(end,2);
            Vss(q,n) = o.V(o.nt).*p.cross_section./1e9;
            V1ss(q,n) = (o.V1(o.nt)+o.V2(o.nt)).*p.cross_section./1e9;
            V2ss(q,n) = (o.V3(o.nt)+o.V4(o.nt)).*p.cross_section./1e9;
        

            %%Linear increase in ELA
            p.da0dt = -7.5/1e4;
            p.tspan=linspace(0,1.5e4,1.5e3+1); %time steps
            p.ic = Rout(end,:);
            dt=p.tspan(2)-p.tspan(1);
            Rout = Rout(end,:);
            for t = p.tspan

                rhs=HMBF_2ODE_2sheets(t,Rout,p);

                R1 = Rout(1) + dt*rhs(1);
                R2 = Rout(2) + dt*rhs(2);

                Rout(1) = R1;
                Rout(2) = R2;
            end


            dVdt = -diff(o.V(2:o.nt))./diff(o.t(2:o.nt)).*p.cross_section./1e9;
            tm = (p.tspan(1:end-1)+p.tspan(2:end))./2;

            MWP_terminate = find(o.R1(1:o.nt)<=0,1);
            MWP_onset = find(o.HS2(1:o.nt)==1,1);

            DVmax(q,n) = max(dVdt(MWP_onset:MWP_terminate));
            DVonset(q,n) = dVdt(MWP_onset);
            DVtime(q,n) = tm(dVdt==DVmax(q,n));
            MWPdur(q,n) = o.t(MWP_terminate) - o.t(MWP_onset);
        end
    end
end

%% Plots
figure(1);set(1,'units','pixels','position',[0 0 1102 562]);
ax1=subplot(1,2,1);
pcolor(h0s,h2s,DVmax./4.6e3);ch=colorbar('southoutside');map = colormap('cool');colormap(ax1,map);shading('flat');caxis([1.5 3]);freezeColors;
hold on;[c1,h1]=contour(h0s,h2s,Vss./1e6,linspace(3e7,4.2e7,5)./1e6,'k','linewidth',4);
clabel(c1,h1,'FontSize',20, 'labelspacing', 300)
hold on;[c3,h3]=contour(h0s,h2s,V2ss./V1ss,[1 2 4],'Color',[0.9 0.9 0.9],'linewidth',4);
clabel(c3,h3,'fontsize',20,'Color',[0.9 0.9 0.9], 'labelspacing', 300)
% hold on;contour(h0s,h2s,Vss,[3.3e7 3.3e7],'r','linewidth',4)
plot(2.4,2.5,'r^','markersize',10,'linewidth',4)
xlabel('$A_1$ ($m^{\frac{1}{2}}$)','fontsize',20,'Interpreter','LaTeX')
ylabel('$A_2$ ($m^{\frac{1}{2}}$)','fontsize',20,'Interpreter','LaTeX')
set(gca,'fontsize',20)
text(0.02,1.0,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)
title('Maximum Loss Rate','fontsize',20)
xlabel(ch,'cm/yr SLE','fontsize',20)


ax2=subplot(1,2,2);
pcolor(h0s,h2s,DVtime);ch2=colorbar('southoutside');map = colormap('winter');colormap(ax2,map);shading('flat');caxis([3e3 5.3e3]);freezeColors;
hold on;[c2,h2]=contour(h0s,h2s,Vss./1e6,linspace(3e7,4.2e7,5)./1e6,'k','linewidth',4);
clabel(c2,h2,'fontsize',20, 'labelspacing', 300)
hold on;[c4,h4]=contour(h0s,h2s,V2ss./V1ss,[1 2 4],'Color',[0.9 0.9 0.9],'linewidth',4);
clabel(c4,h4,'fontsize',20,'Color',[0.9 0.9 0.9], 'labelspacing', 300)
% hold on;contour(h0s,h2s,Vss,[3.3e7 3.3e7],'r','linewidth',4)
plot(2.4,2.5,'r^','markersize',10,'linewidth',4)
xlabel('$A_1$ ($m^{\frac{1}{2}}$)','fontsize',20,'Interpreter','LaTeX')
ylabel('$A_2$ ($m^{\frac{1}{2}}$)','fontsize',20,'Interpreter','LaTeX')
set(gca,'fontsize',20)
text(0.02,1.0,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)
title('Timing','fontsize',20)
xlabel(ch2,'year','fontsize',20)