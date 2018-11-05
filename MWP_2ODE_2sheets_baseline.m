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
p.tspan=[0,4e5]; %time steps

p.plot=0;
options = odeset('RelTol',1e-7,'AbsTol',1e-7,'NonNegative',1);
[time,Rout] = ode23tb(@(t,R) HMBF_2ODE_2sheets(t,R,p),p.tspan,p.ic,options);

R0 = Rout(end,:);
% figure;plot(time./1e3,Rout(:,1)./1e3,time./1e3,Rout(:,2)./1e3)

disp(['Vol is ' num2str(o.V(o.nt).*p.cross_section./1e9 /1e7) ' x 10^7 km^3'])
disp(['Height 1 is ' num2str(p.h0*sqrt(Rout(end,1))) ' m'])
disp(['Height 2 is ' num2str(p.h2*sqrt(Rout(end,2))) ' m'])


%% Deglacial Meltwater pulse with baseline parameters

%Linear increase in ELA
p.da0dt = -7.5/1e4;
p.tspan=linspace(0,1.5e4,1.5e4+1); %time steps
% p.ic = R0;
Rout = R0;
dt=p.tspan(2)-p.tspan(1);
for t = p.tspan

    rhs=HMBF_2ODE_2sheets(t,Rout,p);

    R1 = Rout(1) + dt*rhs(1);
    R2 = Rout(2) + dt*rhs(2);

    Rout(1) = R1;
    Rout(2) = R2;
end

%% Plot Ice6g data

time_start = 18.9;

ncid = netcdf.open('ICE-6G_C_IceThickness_1deg.nc');
thk = netcdf.getVar(ncid,0);
time6g = netcdf.getVar(ncid,1);
lat = netcdf.getVar(ncid,2);
lon = netcdf.getVar(ncid,3);

time6gm = (time6g(1:end-1)+time6g(2:end))./2;

[LAT,LON] = meshgrid(lat,lon);

laur_mask = zeros(size(LAT));
laur_mask(200:305,130:160) = 1;

for n=1:length(time6g)
    lgm_thk = squeeze(thk(:,:,n));
    lgm_thk_laur = lgm_thk.*laur_mask;

    laur_vol(n) = sum(sum(lgm_thk_laur))*(111e3*68e3) / 1e9;
end

subplot(3,1,2)
plot(-time6gm+time_start,diff(laur_vol')./diff(time6g)./1e3./4.6e3,'o','MarkerEdgeColor','w','MarkerFaceColor','k','markersize',10);hold on
box on



%% Plot Tarasov et al., EPSL, 2012 data

load Tarasov_NN5vols_nohead.txt
nn5 = Tarasov_NN5vols_nohead;

nn5(33:2:43,:) = [];

time_nn5_m = (nn5(1:end-1,1)+nn5(2:end,1))./2;
dvdt_nn5 = -diff(nn5(:,2).*25.19)./diff(nn5(:,1))/10;
dvdt_nn5_sig = sqrt(2)*((nn5(1:end-1,3)+nn5(2:end,3))/2)*25.19./(diff(nn5(:,1))*10);

% subplot(3,1,2);plot(time_start + time_nn5_m,dvdt_nn5,'s','MarkerEdgeColor','w','MarkerFaceColor','m','markersize',10);
errorbar(time_start + time_nn5_m,dvdt_nn5,dvdt_nn5_sig,'ms','MarkerEdgeColor','w','MarkerFaceColor','m','markersize',10,'CapSize',10);
xlim([0 10]);hold on;ylim([-1.5 3.5])


%% Calculate some stuff and plot
dVdt = diff(o.V(2:o.nt))./diff(o.t(2:o.nt));
tm = (o.t(2:o.nt-1)+o.t(3:o.nt))./2;
R1 = o.R1(1:o.nt);
R2 = o.R2(1:o.nt);

B1 = o.B1(2:o.nt).*p.cross_section./1e9;B1(R1<=0) = 0;
B2 = o.B2(2:o.nt).*p.cross_section./1e9;B2(R1<=0) = 0;
B3 = o.B3(2:o.nt).*p.cross_section./1e9;B3(R2<=0) = 0;
B4 = o.B4(2:o.nt).*p.cross_section./1e9;B4(R2<=0) = 0;
t = o.t(1:o.nt);

alpha = 1/(1 + (p.h0/p.h2)^2);

ts_plot = [3600 4000 4400 4800 5200];
% clr = parula(length(ts_plot));
q=0;
for t = 2:o.nt
    
    p.a0 = p.a00 + p.da0dt*o.t(t);
    
    R1 = o.R1(t);
    R2 = o.R2(t);
    
    x_i = alpha*((p.h0/p.h2)^2 * R1 - R2 + p.Lc);

    HS = heaviside(R1+R2-p.Lc); %1 when sheets intersect, 0 when separated

    x_i12 = x_i*HS + R1*(1-HS);
    x_i34 = x_i*HS + (p.Lc - R2)*(1-HS);

    x_r1 = (p.h0/(p.s^2))*sqrt(((p.ar-p.a0)/p.beta - p.d0)*p.s + ((p.h0^2) / 4) + R1*(p.s^2)) -...
        ((p.ar-p.a0)/p.beta - p.d0)/p.s - 0.5*(p.h0/p.s)^2;

    x_r2_p = R1 - (1/p.h0^2)*((p.ar-p.a0)/p.beta - p.d0)^2;
    HS2 = heaviside(x_i12-x_r2_p);
    x_r2 = x_r2_p*HS2 + x_i12*(1-HS2);

    x_r3_p = p.Lc - R2 + (1/p.h2^2)*((p.ar-p.a0)/p.beta - p.d0)^2; 
    HS3 = heaviside(x_r3_p-x_i34);
    x_r3 = x_r3_p*HS3 + x_i34*(1-HS3);

    x_r4 = (p.h2/(p.s^2))*sqrt(((p.ar-p.a0)/p.beta - p.d0)*p.s + ((p.h2^2)/4) + R2*(p.s^2)) -...
        ((p.ar-p.a0)/p.beta - p.d0)/p.s - 0.5*(p.h2/p.s)^2;
    
    x_r1 = max([0 x_r1]);
    x_r2 = max([0 x_r2]);
    x_r3 = max([0 x_r3]);
    x_r4 = max([0 x_r4]);
    
    x = linspace(-R1,p.Lc+R2,1e3);
    
    h1 = p.h0.*sqrt(R1-abs(x)); %vialov
    h1(abs(x)>=R1) = 0;

    h2 = p.h2.*sqrt(R2 - abs(x-p.Lc)); %vialov
    h2(p.Lc-x>=R2) = 0;
    h2(p.Lc-x<=-R2) = 0;

    h1(h1<h2) = 0;
    h2(h2<h1) = 0;

    bed(x>0 & x<p.Lc) = p.d0;
    bed(x<=0) = p.d0 + p.s*x(x<=0);
    bed(x>=p.Lc) = p.d0 - p.s*(x(x>=p.Lc)-p.Lc);
    
    H = bed+h1+h2;

    B = p.a0 + p.beta*H;
    B(B>=p.ar) = p.ar;
    
    xp = linspace(-max(o.R1(1:o.nt)),p.Lc+max(o.R2(1:o.nt))+100e3,1e3);
    bedp(xp>0 & xp<p.Lc) = p.d0;
    bedp(xp<=0) = p.d0 + p.s*xp(xp<=0);
    bedp(xp>=p.Lc) = p.d0 - p.s*(xp(xp>=p.Lc)-p.Lc);
    
    intmelt1 = (x>0 & x<R1 & B<p.ar);
    intmelt2 = (x>p.Lc-R2 & x<p.Lc & B<p.ar);
    if(sum(intmelt1) + sum(intmelt2)==0)
        intmelt_area(t) = 0;
        intmelt_avg(t) = p.ar;
    else if(sum(intmelt1)==0)
        intmelt_area(t) = (max(x(intmelt2)) - min(x(intmelt2))).*p.cross_section;
        intmelt_avg(t) = mean(B(intmelt2));
        else
            intmelt_area(t) = (max(x(intmelt1)) - min(x(intmelt1))).*p.cross_section + (max(x(intmelt2)) - min(x(intmelt2))).*p.cross_section;
            intmelt_avg(t) = (mean(B(intmelt1)).*(max(x(intmelt1)) - min(x(intmelt1))).*p.cross_section +...
                mean(B(intmelt2)).*(max(x(intmelt2)) - min(x(intmelt2))).*p.cross_section)./intmelt_area(t);
        end
    end
    
   
    if(sum(o.t(t)==ts_plot)>0)

        q = q+1;   
        figure(1);
        subplot(3,1,1);hold on;
        scatter(x./1e3,H,20,B,'filled');

        cb=colorbar;map = brewermap(32,'RdYlBu');colormap(map);
        caxis([-1 1]);ylabel(cb,'Surface Mass Balance (m/yr)','fontsize',18)
        set(gca,'fontsize',20);xlabel('Distance (km)','fontsize',24);ylabel('Elevation (m)','fontsize',24);
        xlim([-1000 2500]);ylim([-3100 3000]);
    end
end

box on
text(0.01,0.98,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)
area(xp./1e3,bedp,-5000,'FaceColor',[0.690000 0.390000 0.000000],'linewidth',8);hold on
t = o.t(1:o.nt);

text(-840,-1670,[num2str(ts_plot(1)./1000) ' kyr'],'Rotation',-75,'FontSize',20,'Color','k');
text(-700,-1550,[num2str(ts_plot(2)./1000) ' kyr'],'Rotation',-75,'FontSize',20,'Color','k');
text(580,-200,[num2str(ts_plot(3)./1000) ' kyr'],'Rotation',-75,'FontSize',20,'Color','k');
text(380,-200,[num2str(ts_plot(4)./1000) ' kyr'],'Rotation',-75,'FontSize',20,'Color','k');
text(180,-200,[num2str(ts_plot(5)./1000) ' kyr'],'Rotation',-75,'FontSize',20,'Color','k');

figure(1);set(1,'units','pixels','position',[0 0 1002 1202]);
subplot(3,1,2);hold on
plot(t./1e3,-(B1+B4)./4.6e3,'r','linewidth',3);
plot(t./1e3,-(B2+B3)./4.6e3,'b','linewidth',3);
plot(tm(2:end)./1e3,-dVdt(2:end).*p.cross_section./1e9./4.6e3,'k','linewidth',5);hold on
ylabel('Ice loss rate (SLE; cm/yr)','fontsize',20)
set(gca,'fontsize',20);xlabel('Time (kyr)','fontsize',20);
text(0.01,0.98,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)
xlim([0 9]);ylim([-1.5 3.5])


legend('Ice6G Model-Obs','Tarasov et al. Model-Obs','Outer Regions Model','Inner Region Model','Total Model','Location','SouthEast')

subplot(3,1,3)
[ax4,h14,h24] = plotyy(t(2:o.nt)./1e3,intmelt_area(2:o.nt)./1e6./1e6,t(2:o.nt)./1e3,intmelt_avg(2:o.nt));hold on
set(h14,'linewidth',5)
set(h24,'linewidth',5)
xlabel('Time (kyr)','fontsize',20)
ylabel(ax4(1),'Interior area sub-saturation SMB (10^6 km^2)','fontsize',16)
ylabel(ax4(2),'Average interior SMB (m/yr)','fontsize',16)
set(ax4,'fontsize',20,'XLim',[0 9])
set(ax4(1),'YLim',[0 4.5])
set(ax4(2),'YLim',[-6 1.5])
set(gca,'fontsize',16);
text(0.01,0.98,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)
plot(t(find(intmelt_area>0,1)).*ones(1,10)./1e3,linspace(-2,5,10),'-','Color',[0.7 0.7 0.7],'linewidth',3)
plot(t(find(o.HS1(1:o.nt)==0,1)).*ones(1,10)./1e3,linspace(-2,5,10),'--','Color',[0.7 0.7 0.7],'linewidth',3)
plot(t(find(o.R1(1:o.nt)==0,1)).*ones(1,10)./1e3,linspace(-2,5,10),':','Color',[0.7 0.7 0.7],'linewidth',3)


