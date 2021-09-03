% This calculates the error function

function obj = group_model_lss(k,data)

global S1_0 V_0  I1_0 S2_0  I2_0 mu1

time=data.ydata(:,1);   
Aobs= data.ydata(:,2);

E1_0=k(11); U1_0=k(12); E2_0=k(13); U2_0=k(14); 

y0=[S1_0; V_0; E1_0; U1_0; I1_0; S2_0; E2_0; U2_0; I2_0]; 

Amodel = group_predfun(time,k,y0);  
Amodel=mu1*Amodel;
obj = sum((Aobs-Amodel).^2);
