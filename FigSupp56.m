%Symmetric Coevolution
%Jacob Palmer - jacob.palmer@zoo.ox.ac.uk
clc;clear;
close(figure(1));close(figure(2))
%April 19, 2021

%------------------------Leave alone-------------------------------------
n = 3;             %no touchy fishy 
r = zeros(1,n+1);    
r(:,:) = 1;        %Growth rate
tend = 100000;      %Time
Nu1 = 1;
Nu2 = Nu1;
Kn1 = 5;
Kn2 = 5;

%Kmi = 1; %Big number = High MIC.  Low number = low Mic (narrow)
Kmr = 1000;
%Kms = 2; 
Km = 0.05;

DEG = 1;        %0 for no degradation. 1 for degradation

HCE = 1;

%---------------------------Looping------------------------------------
dp = 99;                    %Spectrum steps
dp2 = 199;                   %Steps through abundance
dp3 = 149;                    %Production steps

GamRi0 = 0.0;        %Production rate
GamRiF = 0.45;
GamRS = 0.22198;    %0.22 reported in the text for rounding.  0.22198 actual value ESS calculated and used in this figure. 
GamRr = 0;

cm = 1;           %number of community members 
E = 15;

cy = 0.0001;       %Community abundance T0
Sy = 0.0001;            %Niche competitor abundance T0
Pry = 0;        %Producer Abundance T0 (resident)
Piy0 = 0.000075;
PiyF = 0.000225;

Kmi0 = Km;
KmiF = 3;
Kms = KmiF;

gamS0 = 0;        %Toxin abundance (Niche Competitor) T0
gamr0 = 0;        %Toxin abundance (resident) T0
gami0 = 0;        %Toxin abundance (invader) T0

NO1 = 0.3;        %Niche overlaps
NO2 = 0;

if DEG == 0
    gamD = .75;         %Toxin Degradation
elseif DEG == 1
    gamD = 1;
end
%--------------------------Leave alone-------------------------------------
tspan = [0 tend];
jt = (KmiF-Kmi0)/dp;
gi = (GamRiF-GamRi0)/dp3;
Pc = (PiyF-Piy0)/dp2;
GamIndex = GamRi0:gi:GamRiF;
SigIndex = Kmi0:jt:KmiF;
PyIndex = Piy0:Pc:PiyF;

rc = 0; ic = 0; gc = 0;
endvali = zeros(dp3+1,dp3+1);
GM = zeros(1,dp2+1);
SM = zeros(1,dp2+1);
AM = zeros(1,dp2+1);
Celery = cell(1,dp2+1);
cr = 0;
cS = GamRS;

tic
for Piy = Piy0:Pc:PiyF
    rc = rc + 1
    ic = 0;
    y = [Sy Pry Piy cy gamS0 gamr0 gami0 Nu1 Nu2];
    y0 = y;
%-------------------------Spectrum------------------------------
    for Kmi = Kmi0:jt:KmiF
        ic = ic + 1;
        gc = 0;
            %------------------Production rate focal------------------
            for GamRi = GamRi0:gi:GamRiF
                gc = gc + 1;
                y=y0;
                ci = GamRi; 
                eventfunc = @(t,y) steadystateNutEx(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,DEG,HCE);
                optionsode=odeset('Events',eventfunc, 'NonNegative',1:9);
                if DEG == 1
                    [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
                else
                    [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex_NoDeg(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
                end
                Q = [t,y];
                endvali(gc,ic) = Q(end,4);
            end
    end   
    A = max(max(endvali));
    if A < 1E-9
        A = 0;
        SM(1,rc) = NaN;
        GM(1,rc) = 0;
    else
        [GamMax,SigMax] = find(endvali==A);
        SigMax = SigIndex(SigMax);
        GM(1,rc) = GamIndex(GamMax);
        SM(1,rc) = (1 - SigMax/KmiF) / (1 - Kmi0/KmiF);
    end
    AM(1,rc) = A;
    Celery{1,rc} = endvali;
end

save('C:\Users\jdpal\Desktop\Local\Data\FigSupp56_075-225_GamRiF045_dp3-149_dp2-199','Celery')

figure(1)
yyaxis left
plot(PyIndex,SM)
xlim([0.000075 0.000225])
xlabel('Starting Abundance')
xticks([0.0001 0.00015 0.0002])
ylabel('\sigma_{AC}')
ylim([-0.05 1.05])
hold on
yyaxis right
plot (PyIndex,GM)
ylabel('\gamma_{A}')
ylim([0 0.45])
hold off

% figure(2)
% plot(PyIndex,AM)
% ylim([0 max(max(AM))])
% ylabel('Final Abundance - A')
% xlabel('Starting Abundance - A')
    
tim = toc;
persec = (dp*dp2*dp3) / tim
% load splat
% sound(y,Fs)
