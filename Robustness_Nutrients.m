%Symmetric Coevolution
%Jacob Palmer - jacob.palmer@zoo.ox.ac.uk
clc;clear;
%April 19, 2021

%------------------------Leave alone-------------------------------------
n = 3;             %no touchy fishy 
r = zeros(1,n+1);    
r(:,:) = 1;        %Growth rate
tend = 100000;      %Time
Nu1 = 0.5;
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
dp3 = 149;                    %Production steps

GamRi0 = 0;        %Production rate
GamRiF = 0.95;
GamRs0 = GamRi0;
GamRsF = GamRiF;
GamRr = 0;

cm = 1;           %number of community members 
E = 15;

cy = 0.0001;       %Community abundance T0
Sy = 0.0001;            %Niche competitor abundance T0
Pry = 0;        %Producer Abundance T0 (resident)
Piy = Sy;

Kms0 = Km;
KmsF = 3;
Kmi0 = Km;
KmiF = KmsF;

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
gr = (GamRsF-GamRs0)/dp3;
GamIndex = GamRi0:gi:GamRiF;

rc = 0; ic = 0; gc = 0; gt = 0;
endvali = zeros(dp3+1,dp3+1);
endvalS = zeros(dp3+1,dp3+1);
SigmaDiag = zeros(1,dp+1);
ESSr = zeros(1,dp+1);
y = [Sy Pry Piy cy gamS0 gamr0 gami0 Nu1 Nu2];
y0 = y;
cr = 0;
eCelery = cell(1,dp+1);

tic
%-------------------------Spectrum------------------------------
for Kmi = Kmi0:jt:KmiF
    Kms = Kmi;
    ic = ic + 1
    gc = 0;
        %------------------Production rate focal------------------
        for GamRi = GamRi0:gi:GamRiF
            gc = gc + 1;
            gt = 0;
            %---------------Production rate resident---------------
            for GamRS = GamRs0:gr:GamRsF
                gt = gt + 1;
                y=y0;
                ci = GamRi; 
                cS = GamRS;
                eventfunc = @(t,y) steadystateNutEx(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,DEG,HCE);
                optionsode=odeset('Events',eventfunc, 'NonNegative', 1:9);
                if DEG == 1
                    [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
                else
                    [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex_NoDeg(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
                end
                Q = [t,y];
                endvali(gc,gt) = Q(end,4);
            end
        end
        eCelery{1,ic} = endvali;
end   

save('C:\Users\jdpal\Desktop\Local\Data\FigSupplement2_Nutrients_Nu05','eCelery')

tim = toc
persec = (dp*dp3*dp3) / tim
% load splat
% sound(y,Fs)
