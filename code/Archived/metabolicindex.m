load(['trawl_data_2008.mat']);
%O2 from trawl data is in ml/l - may need to be converted to umol/kg

gas_const = 8.31;
partial_molar_vol = 0.000032;
kelvin=273.15;
Boltz= 0.000086173324;


%calculate percent saturation for O2 - assumes  units of mL O2/L.
% Input:       S = Salinity (pss-78)
%              T = Temp (deg C) ! use potential temp
%depth is in meters
%[umole/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)


[SA, in_ocean] = gsw_SA_from_SP(sal,depth,lon,lat); %absolute salinity for pot T calc

pt = gsw_pt_from_t(SA,temp,depth); %potential temp at a particular depth

CT = gsw_CT_from_t(SA,temp,depth); %conservative temp

sigma0 = gsw_sigma0(SA,CT);


o2_umolkg = o2.*44660./(sigma0+1000);

[O2_Sat0] = gsw_O2sol_SP_pt(sal,pt); %relies on o2 in umol/kg
%= o2satv2a(sal,pt); %uses practical salinity and potential temp - solubity at p =1 atm

press = exp(depth.*10000.*partial_molar_vol./gas_const./(temp+kelvin));

O2_satdepth=O2_Sat0.*press;

%solubility at p=0
sol0=O2_Sat0./0.209;

sol_Dep = sol0.*press;

po2=o2_umolkg./sol_Dep;

%met index phi

%assuming its a cod 
%Ao = 3.11E-14;
%Eo =  0.8736;
%B=500;% size 
%N=-0.208; %assuming its a cod 

%assuming its a giant monster sablefish 
Ao = 1.16625e-13;
Eo =  0.8736;
B= 10000;% size 
N=-0.208; %assuming its a cod 

phi = B^N.*Ao.* po2 ./ exp(-1*Eo./ (Boltz * (temp+kelvin)));


save(['trawl_data_2008_phi.mat']);

