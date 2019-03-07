% this function is used to calculate the Pressure and temperature related
% parameters in the superheat region 
% Tsat = saturate temperature
% h = heat transfer coeficient 
% Cp = Specific heat
% v = volumn?
function [Tsat,h,Cp,v]=hsh(P,T)
T=T-273;
% [Tsat]=hsat(P);
[Tsat,dum1,dum2]=hsat(P);
Tsat=Tsat-273;
if P>2.955 && P<6.522
    hgg=-1.347244*(P-2.999)^2.0-2.326913*(P-2.999)+2803.35;
elseif P>6.522 && P< 16.497
    hgg=-0.9219176*(P-9.0)^2.0-16.38835*(P-9.0)+2742.03;
end;
h=hgg+(4.5*P/(7.4529E-6*T^3-P^2)^0.5+0.28*exp(-0.008*(T-162))-100/T+2.225)*(T-Tsat);
h=h*1000;
if P>3.932 && P< 8.996
    vgg=(2.868721*P^1.252148+3.8)^(-1);
elseif P>8.996 && P< 14.628
    vgg=(0.5497653*P^1.831182+18.111)^(-1);
elseif P>14.628 && P< 18.210
    vgg=(8.5791582E-3*P^3.176484+50.0)^(-1);
end;
v=vgg+(0.000466/P-(0.12/(T+100)-0.00106)*P^0.1/(1.96E-8*(T+8)^4-P^2)^0.5)*(T-Tsat);
if P>2.391 && P< 5.661
Cpg=0.3187082*P^1.110271+2.3;
elseif P>5.661 && P< 9.458
    Cpg=0.064275995*P^1.766106+3.12;
elseif P>9.458 && P<12.9
    Cpg=3.8011048E-3*P^2.816897+4.40;
elseif P>12.9 && P<16.309
    Cpg=0.1876175*exp(0.2466925*P)+5.0;
elseif P>16.309 && P<18.743
    Cpg=7.620756E-3*exp(0.4117289*P)+9.2;
end;
Cp=Cpg-(0.011*P/(0.00014*(T+8)^2-P)^1.5+1.5E-8*(655-T)^2.1)*P*(T-Tsat);
Cp=Cp*1000;
Tsat=Tsat+273;
return;   

    
    