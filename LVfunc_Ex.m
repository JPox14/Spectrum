function dydt = LVfunc_Ex(t, y, r, cS, ci, cr, gamD,Kn1,Kn2,Kmi,Kmr,Kms,Km,NO1,NO2,E,HCE)

dydt = zeros(1,9);
dydt = dydt';

if y(1) < 1E-9 
    y(1) = 0;
end

if y(2) < 1E-9
    y(2) = 0;
end

if y(3) < 1E-9
    y(3) = 0;
end

if y(4) < 1E-9
    y(4) = 0;
end

if y(5) < 1E-9
    y(5) = 0;
end

if y(6) < 1E-9
    y(6) = 0;
end

if y(7) < 1E-9
    y(7) = 0;
end

if y(8) < 1E-9
    y(8) = 0;
end

if y(9) < 1E-9
    y(9) = 0;
end

%Monod
y7 = y(1)*gamD*(y(7)^HCE/(y(7)^HCE+Km^HCE));
y7I = y(4)*gamD*(y(7)^HCE/(y(7)^HCE+Kmi^HCE));

y6 = y(1)*gamD*(y(6)^HCE/(y(6)^HCE+Km^HCE));
y6r = y(4)*gamD*(y(6)^HCE/(y(6)^HCE+Kmr^HCE));

y52 = y(2)*gamD*(y(5)^HCE/(y(5)^HCE+Km^HCE));
y53 = y(3)*gamD*(y(5)^HCE/(y(5)^HCE+Km^HCE));
y5s = y(4)*gamD*(y(5)^HCE/(y(5)^HCE+Kms^HCE));

% y7 = y(1)*(1/(1+(Km/y(7))^HCE));
% y7I = y(4)*(1/(1+(Km/y(7))^HCE));
% 
% y6 = y(1)*(1/(1+(Km/y(6))^HCE));
% y6r = y(4)*(1/(1+(Km/y(6))^HCE));
% 
% y52 = y(2)*(1/(1+(Km/y(5))^HCE));
% y53 = y(3)*(1/(1+(Km/y(5))^HCE));
% y5s = y(4)*(1/(1+(Km/y(5))^HCE));

%Morrison
% y7 = (((y(1)+y(7)+Km) - sqrt( (y(1)+y(7)+Km)^2 - 4*(y(1)*y(7)) ))  /  2*y(1));
% y7I = (((y(4)+y(7)+Kmi) - sqrt( (y(4)+y(7)+Kmi)^2 - 4*(y(4)*y(7)) ))  /  2*y(4));
% 
% y6 = (((y(1)+y(6)+Km) - sqrt( (y(1)+y(6)+Km)^2 - 4*(y(1)*y(6)) ))  /  2*y(1));
% y6r = (((y(4)+y(6)+Kmr) - sqrt( (y(4)+y(6)+Kmr)^2 - 4*(y(4)*y(6)) ))  /  2*y(4));
% 
% y52 = (((y(2)+y(5)+Km) - sqrt( (y(2)+y(5)+Km)^2 - 4*(y(2)*y(5)) ))  /  2*y(2));
% y53 = (((y(3)+y(5)+Km) - sqrt( (y(3)+y(5)+Km)^2 - 4*(y(3)*y(5)) ))  /  2*y(3));
% y5s = (((y(4)+y(5)+Kms) - sqrt( (y(4)+y(5)+Kms)^2 - 4*(y(4)*y(5)) ))  /  2*y(4));


%Invader Toxin
    dydt(7) =ci*y(3)*(y(8)/(y(8)+Kn1)) - ((y7 + y7I));
%Resident Toxin
    dydt(6) = cr*y(2)*(y(8)/(y(8)+Kn1)) - ((y6 + y6r));
%Niche Competitor Toxin
    dydt(5) = cS*y(1)*(y(8)/(y(8)+Kn1)) - ((y52 + y53 + y5s));

%Invader
dydt(3)=y(3)*((y(8)/(y(8)+Kn1))*r(3)*(1-ci) - E*(y(5)^HCE/(y(5)^HCE+Km^HCE)));

%Resident
dydt(2)=y(2)*((y(8)/(y(8)+Kn1))*r(2)*(1-cr) - E*(y(5)^HCE/(y(5)^HCE+Km^HCE)));

%Niche Competitor
dydt(1)=y(1)*((y(8)/(y(8)+Kn1))*r(1)*(1-cS) - (E*(y(6)^HCE/(y(6)^HCE+Km^HCE)) + E*(y(7)^HCE/(y(7)^HCE+Km^HCE))));

%Community
dydt(4)=y(4)*((y(9)/(y(9)+Kn2))*r(4) + (y(8)/(y(8)+Kn1))*r(4)*NO1 - E*(y(5)^HCE/(y(5)^HCE+Kms^HCE)) - E*(y(6)^HCE/(y(6)^HCE+Kmr^HCE)) - E*(y(7)^HCE/(y(7)^HCE+Kmi^HCE)));

%Nutrients
dydt(8)=-(y(8)/(y(8)+Kn1))*(y(1)+y(2)+y(3)+y(4)*NO1);
dydt(9)=-(y(9)/(y(9)+Kn2))*(y(4)+y(3)*NO2+y(2)*NO2+y(1)*NO2);

%Nurients released back into environment
% dydt(8) = dydt(8) + (y(8)/(y(8)+Kn1))*(y(7)*E*(y(1) + SPi*y(4)*NO1)) + (y(6)*E*(y(1) + SPr*y(4)*NO1)) + (y(5)*E*(y(2) + y(3) + SS*y(4)*NO1))*.8;
% dydt(9) = dydt(9) + (y(9)/(y(9)+Kn2))*(y(5)*E*(SS*y(4))) + (y(6)*E*(SPr*y(4))) + (y(7)*E*(SPi*y(4)))*.8;

end

