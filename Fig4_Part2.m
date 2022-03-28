%eCelery = load(filepath);
%eCelery = eCelery.eCelery;
% dp = 99;
% dp3 = 299;
% GamRi0 = 0;        
% GamRiF = 0.4;
% gi = (GamRiF-GamRi0)/dp3;
% GamIndex = GamRi0:gi:GamRiF;

xi = 0;
total = 0;
SigmaDiag1 = zeros(1,dp+1);
ESSr1 = zeros(1,dp+1);
toad = 0;


for x = 1:1:dp+1
    endvali = eCelery{1,x};
    xi = xi + 1;
    for k = 1:1:dp3+1
        for m = k:1:dp3
            %k=1
            if (k == 1)
                if (endvali(k,k) > endvali(k+m,k))
                    total = total + 1;
                end
            end
            if (k > 1) && (endvali(k,k) > endvali(m+1,k))
                total = total + 1;
            end
        end
        
        if k > 1
            for M = k:-1:2
                if (endvali(k,k) > endvali(M-1,k))
                    total = total + 1;
                end
            end
        end
        
        if total == 0
            toad = 0;
        elseif toad < total
            toad = total;
            SigmaDiag1(1,xi) = GamIndex(k);
        elseif (toad == total) && (k > 1)
            SigmaDiag1(1,xi) = (GamIndex(k) + GamIndex(k-1)) / 2;
            toad = total;
        end
    total = 0;
    end
    toad=0;
end

% for y = dp:-1:2
%     if (SigmaDiag1(1,y) == 0)
%         SigmaDiag1(1,y) = (SigmaDiag1(1,y-1) + SigmaDiag1(1,y+1)) / 2;
%     end
% end

figure(6)
plot(flip(SigmaDiag1))
xlabel('\sigma_{AC} , \sigma_{BC}')
ylabel('\gamma_{ESS}')
xticks([0 dp*0.25 dp*0.5 dp*0.75 dp+1*1])
xticklabels({'0','','','','1'})
xlim([1 dp+2])
ylim([0 (max(SigmaDiag1)+.01)])
% 