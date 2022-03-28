function [value,isterminal,direction] = steadystateNutEx(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,DEG,HCE)


%must decide what to do below
if DEG == 1
    dy = LVfunc_Ex(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE);
else
    dy = LVfunc_Ex_NoDeg(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE);
end
    
%Record when steady state is reached
SS = [abs(dy(8)), abs(dy(9))];
SS1 = max(SS);
val1 = SS1 - 0.00000005;

value = val1;
isterminal = 1;
direction = -1;

end