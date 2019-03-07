% this function is used to calculate the parameters under saturation
% condition. the subscripting symbol f stands for fluid, and g stands for
% gas. 
% p = inlet pressure
% Tsat = saturation temperature
% h = heat transfer coefficient
% k = Thermal conductivity (W/m C)
% mu = viscocity ()
% Pr = Prantl number

function [Tsat,hf,hg,kf,kg]=hsat(P)

if P>1.676 && P<8.511
    Tsat=185.0779*P^0.2323217-5.0;
elseif P>8.511 && P<17.690
    Tsat=195.1819*P^0.2241729-16.0;
end;
Tsat=Tsat+273;

if P>4.0200 && P<9.964
    a=373.7665; b=0.4235532; c=415.00;
elseif P>9.964  && P<16.673
    a=75.38673; b=0.8282384; c=900;
end;
hf=a*P^b+c;
hf=hf*1000;

if P>2.955 && P<6.522
    a=-1.347244; b=-2.999; c=-2.326913; d=2803.35;
elseif P>6.5222  && P<16.497
    a=-0.9219176; b=-9.0; c=-16.38835; d=2742.03;
end;
hg=a*(P+b)^2+c*(P+b)+d;
hg=hg*1000;

PP0=[3.9776,4.6941,5.5052,6.4191,7.4449];
kkf=[0.616,0.603,0.589,0.574,0.558];
kkg=[49.5,52.8,56.6,60.9,66.0]*1.0E-3;
Prff=[0.859,0.866,0.882,0.902,0.932];
Prgg=[1.39,1.43,1.48,1.54,1.61];

kf=interp1(PP0,kkf,P);
kg=interp1(PP0,kkg,P);
% Prf=interp1(PP0,Prff,P);
% Prg=interp1(PP0,Prgg,P);
% 
% 
% if P>3.948 && P<9.514
%     muf=141.5415-25.91353*log(P);
%     muf=muf*1.0E-6;
% end;
% 
% if P>2.207 && P<5.480
%     mug=(3.375163*P^(0.3916208)+11.8)*1.0E-6;
% elseif P>5.480 && P<9.585
%     mug=(0.9169410*P^(0.7644731)+15.0)*1.0E-6;
% end;



return;


    

    
    