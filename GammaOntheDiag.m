y =95;
x = dp+2 - y;
 %x=2;
mat4 = eCelery{1,x};
wr1 = diag(mat4)';
I1 = mat4./wr1;

% Create a matrix where each row is a color triple
hexvals1 = ['#3B4244';'#FFFFFF';'#EA696D']; %black - white - red
V1 = hex2rgb(hexvals1);
X1 = [0 0.5 1];
Xq = linspace(0, 1, 255);
Vq1 = interp1(X1, V1, Xq, 'pchip');

figure(10)
imagesc(I1)
caxis([0.9999999 1.0000001])
colormap(Vq1)
cb = colorbar;
cb.Label.String = '(No Invasion              Invasion)';
set(gca, 'YDir','normal','XDir','normal')
xlabel('\gamma_{B}')
ylabel('\gamma_{A}')
xticks([dp3*0.01 dp3*0.25 dp3*0.525 dp3*0.75 dp3*1])
xticklabels({'0%','','50%','','95%'})
yticks([0 dp3*0.25 dp3*0.525 dp3*0.75 dp3*1])
yticklabels({'0%','','50%','','95%'})