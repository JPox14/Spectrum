function dydt = LVfunc_Ex_NoDeg(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE)

dydt = zeros(1,9);
dydt = dydt';

%Invader Toxin
    dydt(7) =ci*y(3)*(y(8)/(y(8)+Kn1)) - gamD*y(7);
%Resident Toxin
    dydt(6) = cr*y(2)*(y(8)/(y(8)+Kn1)) - gamD*y(6);
%Niche Competitor Toxin
    dydt(5) = cS*y(1)*(y(8)/(y(8)+Kn1)) - gamD*y(5);

%Invader
dydt(3)=y(3)*((y(8)/(y(8)+Kn1))*r(3)*(1-ci) - E*(y(5)^HCE/(y(5)^HCE+Km^HCE)));

%Resident
dydt(2)=y(2)*((y(8)/(y(8)+Kn1))*r(2)*(1-cr) - E*(y(5)^HCE/(y(5)^HCE+Km^HCE)));

%Niche Competitor
dydt(1)=y(1)*((y(8)/(y(8)+Kn1))*r(1)*(1-cS) - (E*(y(6)^HCE/(y(6)^HCE+Km^HCE)) + E*(y(7)^HCE/(y(7)^HCE+Km^HCE))));

%Community
% dydt(4)=y(4)*((y(9)/(y(9)+Kn2))*r(4)  + (y(8)/(y(8)+Kn1))*r(4)*NO1 - E*(y(5)^HCE/(y(5)^HCE+Kms^HCE)) - E*(y(6)^HCE/(y(6)^HCE+Kmr^HCE)) - E*(y(7)^HCE/(y(7)^HCE+Kmi^HCE)));

if y(4) < 0.00001
    %y(4) = 0;
    dydt(4) = -10;
else
    dydt(4)=y(4)*((y(9)/(y(9)+Kn2))*r(4)  + (y(8)/(y(8)+Kn1))*r(4)*NO1 - E*(y(5)^HCE/(y(5)^HCE+Kms^HCE)) - E*(y(6)^HCE/(y(6)^HCE+Kmr^HCE)) - E*(y(7)^HCE/(y(7)^HCE+Kmi^HCE)));
end
    
%Nutrients
dydt(8)=-(y(8)/(y(8)+Kn1))*(y(1)+y(2)+y(3)+y(4)*NO1);
dydt(9)=-(y(9)/(y(9)+Kn2))*(y(4)+y(3)*NO2+y(2)*NO2+y(1)*NO2);

%Nutrients released back into environment
%dydt(8) = dydt(8) + (y(7)*(y(1) + SPi*cm*y(4)*NO1)) + (y(6)*(y(1) + SPr*cm*y(4)*NO1)) + (y(5)*(y(2) + y(3) + SS*cm*y(4)*NO1));
%dydt(9) = dydt(9) + y(5)*(SS*cm*y(4)) + y(6)*(SPr*cm*y(4)) + y(7)*(SPi*cm*y(4));


end