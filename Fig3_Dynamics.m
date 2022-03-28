%Lotka-Volterra 3 species Aug26 2020
%Jacob Palmer - jacob.palmer@zoo.ox.ac.uk
%clc;clear;
%close(figure());
%Lotka_Volterra_Ex_Dynamics.m

%Play with GamR to see how a local success can result in less biomass but
%total removal of the Niche Competitor.  

%When E = 100, then toxin production can result in greater biomass than no
%toxin production, while also removing the NC.  Otherwise, its always
%better to not make a toxin, based on total biomass. 


%------------------------Leave alone-------------------------------------
n = 3;             %no touchy fishy
r = zeros(1,n+1);
r(:,:) = 1.0;        %Growth rate
Nu1 = 1;
Nu2 = Nu1;
Kn1 = 5;
Kn2 = 5;

%Toxin-receptor affinity  -- Big number (3) = low affinity.  Low number (0.05) = High affinity
Kmi = 0.0795;  %focal strain -> community
Kmr = 0.05; %ignore
Kms = 0.0795; %conspecific strain -> community
Km = 0.05;  %Fixed High affinity between focal strain and conspecific

%Km = 0.0798 : GamR-ESS = 0.16577 (Sigma = 0.98989, rounded to 0.99 in text
%Km = 3.0 :    GamR-ESS = 0.22315 (Sigma = 0).


HCE = 1;        %Hill coefficient.  

DEG = 1;        %0 for no degradation. 1 for degradation

%-------------------------Factors to play with----------------------------

GamR = 0.16577;        %production rate

GamRi = GamR;
GamRr = 0;        %ignore
GamRS = GamR;

cm = 1;           %number of community members
E = 15;          %Killing Efficiency of the toxin

cy = 0.0001;         %Community abundance T0
Sy = cy;             %Niche competitor abundance T0
Pry = 0;             %Producer Abundance T0 (resident)
Piy = 0.0001;        %Producer Abundance T0 (invader)

gamS0 = 0;        %Toxin abundance (Niche Competitor) T0
gamr0 = 0;        %Toxin abundance (resident) T0
gami0 = 0;        %Toxin abundance (invader) T0

NO1 = 0.3;        %Niche overlap
NO2 = 0;

%------------------------Factors to leave alone (for now)------------------

if DEG == 0
    gamD = .75;         %k -- Toxin Degradation
elseif DEG == 1
    gamD = 1;           %Theta -- Absorption term
end

ci = GamRi; 
cr = GamRr; 
cS = GamRS;

tend = 100000;      %Time

%--------------------------Leave alone-------------------------------------
y = [Sy Pry Piy cy gamS0 gamr0 gami0 Nu1 Nu2];
y0 = y;
tspan = [0 tend];

%--------------------------------------------------------------------------

tic
%-------------------------------------Numerical solution-----------------------------
eventfunc = @(t,y) steadystateNutEx(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,DEG,HCE);
optionsode=odeset('Events',eventfunc,'NonNegative',1:9);
if DEG == 1
    [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
else
    [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex_NoDeg(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
end
M = [t,y];
b = M(end,4);
d = M(end,3);
e = M(end,5);
%b = b+d;
N = M(:,[2:5]);
N1 = M(:,2);
N2 = M(:,3);
N3 = M(:,4);
N4 = M(:,5);
O = M(:,9);
o = M(:,10);
Toxs = M(:,6);
Toxr = M(:,7);
Toxi = M(:,8);
O = O./Nu1;
o = o./Nu2;

cmax = max(N4);
fmax = max(N1) + max(N2) + max(N3);
N1 = N1/fmax;
N2 = N2/fmax;
N3 = N3/fmax;
N4 = N4/cmax;

toxmax = max(Toxs) + max(Toxr) + max(Toxi);
Toxs = Toxs/toxmax;
Toxr = Toxr/toxmax;
Toxi = Toxi/toxmax;

b = num2str(b);

if DEG == 1
    figure(8)
else
    figure(9)
end

plot(t, N1, 'color',[0, 0.4470, 0.7410],'LineWidth',2)
title(b)
hold on
plot(t,N3,'color',[0.8588, 0.2667, 0.2157],'LineWidth',2)
plot(t,Toxs,'--','color',[0, 0.4470, 0.7410],'LineWidth',1)
plot(t,Toxi,'--','color',[0.8588, 0.2667, 0.2157],'LineWidth',1)
plot(t,N4,'color',[0.9, 0.9, 0.0],'LineStyle','-','LineWidth',2)
plot(t,O,'color',[0.99, 0.5250, 0.0],'LineStyle','-','LineWidth',1)
plot(t,o,'color',[0.05 0.6157 0.3451],'LineStyle','-','LineWidth',1)

%legend('Niche Competitor','Focal species','Community','Toxin','Nutrient 1','Nutrient 2','Location','northwest','numcolumns',n)
xlabel('Time (Hours)')
ylabel('Relative Abundance')
ylim([0 1.1]);
if DEG == 1
    xlim([0 110]);
else
    xlim([10 100]);
end
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})


hold off

toc