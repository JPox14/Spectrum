
endvalI = zeros(dp+1,dp+1);
endvalB = zeros(dp+1,dp+1);
sc = 0;

%Use from ReMakeSigmaDiag
SigmaDiag = SigmaDiag1;

tic
for Kms = Kmi0:jt:KmiF
sc = sc+1
ic = 0;
    for Kmi = Kmi0:jt:KmiF
        ic = ic+1;
        ci = SigmaDiag(1,sc);
        cS = SigmaDiag(1,sc);
        y=y0;
        eventfunc = @(t,y) steadystateNutEx(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,DEG,HCE);
        optionsode = odeset('Events',eventfunc, 'NonNegative',1:9);
        if DEG == 1
            [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
        else
            [t,y,te,ye,ie] = ode45(@(t,y) LVfunc_Ex_NoDeg(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE),tspan, y0,optionsode);
        end
        Q = [t,y];
        endvalI(ic,sc) = Q(end,4);
        endvalB(ic,sc) = Q(end,2);
        %build the logic to find the max value in each spectrum column
    end
end

%save('C:\Users\jdpal\Desktop\Local\Fig2PIP_Part3_eCelery_E15_GamRiF-Tight_Ext','eCelery')

iPRO = endvalI./diag(endvalI)';

% Create a matrix where each row is a color triple
hexvals1 = ['#3B4244';'#FFFFFF';'#EA696D']; %black - white - red
hexvals2 = ['#FFFFFF';'#FCF3F1';'#E11C1C']; %white - red
hexvals3 = ['#FFFFFF';'#F4E1FF';'#401858']; %white - purple
V1 = hex2rgb(hexvals1);
V2 = hex2rgb(hexvals2);
V3 = hex2rgb(hexvals3);

% And now make a vector of what range each color should be at (i.e. this vector 
% defines the spacing of the colours, they need not be regularly/equally spaced):
X1 = [0 0.5 1];
X2 = [0 0.1 1];
X3 = [0 0.1 1];
% 
% And finally you can create the entire map with one simple interpolation:
Xq = linspace(0, 1, 255);
Vq1 = interp1(X1, V1, Xq, 'pchip');
Vq2 = interp1(X2, V2, Xq, 'pchip');
Vq3 = interp1(X3, V3, Xq, 'pchip');

figure(5)
imagesc(iPRO)
caxis([0.999995 1.000005])
colormap(Vq1)
cb = colorbar;
cb.Label.String = '(No Invasion              Invasion)';
set(gca, 'YDir','reverse','XDir','reverse')
xlabel('\sigma_{BC}')
ylabel('\sigma_{AC}')
xticks([0 dp*0.25 dp*0.5 dp*0.75 dp*1])
xticklabels({'1','','','',''})
yticks([0 dp*0.25 dp*0.5 dp*0.75 dp*1])
yticklabels({'1','','','','0'})

 tim = toc;
 persec1 = (dp*dp) / tim
load splat
sound(y,Fs)