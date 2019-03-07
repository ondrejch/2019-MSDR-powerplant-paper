% This function is used to calculate the specific enthalpy in
% sub-cooling region.
% Inputs:
%   P = the inlet pressure
%   T =  inlet temperature
% Outputs:
%   Ps  =  saturation pressure
%   h   =  specific enthalpy 
%   Cp  =  specific heat
%   rho =  steam density
function [Ps,h,Cp,rho]=hsub(P,T);

T=T-273;
if T>89.965 & T<179.781
    Ps=((T+57.0)/236.2315)^5.602972;
elseif T>139.781 & T<203.662
    Ps=((T+28.0)/207.9248)^4.778504;
elseif T>203.662 & T<299.407
    Ps=((T+5.0)/185.0779)^4.304376;
elseif T>288.407 & T<355.636
    Ps=((T+16.0)/195.1819)^4.460843;
elseif T>355.636 & T<373.253
    Ps=((T+50.0)/227.2963)^4.960785;
end;
if Ps>0.075 & Ps<0.942
    h=912.1779*Ps^0.2061637-150.0;
elseif Ps>0.942 & Ps<4.02
       h=638.0621*Ps^0.2963192+125.0;
elseif Ps>4.020 & Ps<9.964
h=373.7665*Ps^0.4235532+415.0;
elseif Ps>9.964 & Ps<16.673
h=75.38673*Ps^0.8282384+900.0;
end;
h=h*1000;
if Ps>0.03 & Ps<0.671
    cpf=0.247763*Ps^0.5704026+4.15;
elseif Ps>0.671 & Ps<2.606
    cpf=0.1795305*Ps^0.8967323+4.223;
elseif Ps>2.606 & Ps<6.489
    cpf=0.0935984*Ps^1.239114+4.340;
elseif Ps>6.489 & Ps< 11.009
    cpf=0.01068888*Ps^2.11376+4.740;
elseif Ps>11.009 & Ps< 14.946
    cpf=1.0E-4*1.333058*Ps^3.707294+5.480;
elseif Ps>14.946 & Ps< 18.079
    cpf=1.0E-3*6.635658*(Ps-10.0)^3.223323+7.350;
end;
Cp=cpf+(0.0018-76/((364-T)^1.8))*(P-Ps);
Cp=Cp*1000;

if Ps>0.075 & Ps<1.000
    rhof=(1.2746977*1.0E-4*Ps^0.4644339+0.001)^(-1);
elseif Ps>1.00 & Ps<3.88
    rhof=(1.0476071*1.0E-4*Ps^0.5651090+0.001022)^(-1);
elseif Ps>3.880 & Ps<8.840
    rhof=(3.2836717*1.0E-5*Ps+1.12174735*1.0E-3)^(-1);
elseif Ps>8.840 & Ps< 14.463
    rhof=(3.3551046*1.0E-4*exp(5.8403566*1.0E-2*Ps)+0.00085)^(-1);
elseif Ps>14.463 & Ps< 18.052
    rhof=(3.1014626*1.0E-8*Ps^3.284754+0.001430)^(-1);
end;
rho=rhof+(170/(375-T)-0.2)*(P-Ps);
return;


    

    
    