% This is the main file that is used to obtain the mcmc estimate. 
clc
clear model data params options
close all

load cdata.mat
data.ydata=groupdatacumulative;

global S1_0 V_0  I1_0 S2_0 I2_0 mu1 mu2

mu1=1/3; mu2=mu1;
%===============
% Population data:
% Total population in ontaio at 2109: 14,570,000
% Based on 2016 census data, the distribution of age group between 0-14 is 16.4%
% In 2016, the population between 0-14 is  2,207,970; the population
% between 15-19 is 811,670
% We can estimate the distribution of age group 15-19 is 6.03%
% The percentage of age group between 0-19 is: (16.4+6.03)%
%================
young_p=(16.4+6.03)/100;
old_p=1-young_p; 
young_ip=407/(407+295+277+222+230+157+89+53+27);
old_ip=1-young_ip;
cumulative_n=184635;

S1_0= 14570000*old_p-cumulative_n*old_ip; % Initial population of susceptible individuals above age 20 at Jan. 1, 2021
V_0= S1_0*1e-5;  
I1_0= 2476*old_ip; 
S2_0= 14570000*young_p-cumulative_n*young_ip; % Initial population of susceptible individuals below age 20 at Jan. 1, 2021
I2_0= 2476*young_ip;

format shortEng
format compact

opts = optimset('fminsearch');
opts.TolX = 1.e-6;
opts.MaxFunEvals = 1e5;
opt.MaxIter=1e5;

para_in=[1e-8; 1e-8; 1e5; 1e-7; 0.95; 1e-7; 1e-5; 1e-8; 1e-8; 1e-5; 500; 100; 500; 100]; 

fun_k=@(k) group_model_lss(k,data); 
LB=[0; 0; 0; 0; 0.5; 0; 0; 0; 0; 0; 0; 0; 0; 0]; 
UB=[1e-6; 1e-6; 1; 1e-6; 1; 1e-6; 1e-3; 1e-6; 1e-6; 1e-3; 3000; 2000; 3000; 2000];
[k0,ss0] = fminsearchbnd(fun_k,para_in,LB,UB,opts)

mse = ss0/(length(data.ydata(:,1))-length(para_in)); 

params = {
    {'beta11',  k0(1), 0,1e-6}
    {'beta12',  k0(2), 0,1e-6}
    {'omega',   k0(3), 0,1}
    {'beta11v', k0(4), 0,1e-6}
    {'epsilon', k0(5), 0.5,1}
    {'beta12v', k0(6), 0,1e-6}
    {'delta1',  k0(7),0,1e-3}
    {'beta21',  k0(8),0,1e-6}
    {'beta22',  k0(9),0,1e-6}
    {'delta2',  k0(10),0,1e-3}
    {'E1_0',    k0(11),0,3000}
    {'U1_0',    k0(12),0,2000}
    {'E2_0',    k0(13),0,3000}
    {'U2_0',    k0(14),0,2000}
    };


model.ssfun = @group_model_lss;
model.sigma2 = mse;


options.nsimu = 1000;
options.burnintime=700;
options.updatesigma = 1;
options.method='mh';

[results,chain,s2chain] = mcmcrun(model,data,params,options);

figure(1); 
mcmcplot(chain,[],results,'chainpanel')

figure(2); 
mcmcplot(sqrt(s2chain),[],[],'dens',2)
title('error std')

chainstats(chain,results)

para_m=mean(chain);     
E1_0m=para_m(11); U1_0m=para_m(12); E2_0m=para_m(13); U2_0m=para_m(14); 

y0_m=[S1_0; V_0; E1_0m; U1_0m; I1_0; S2_0; E2_0m; U2_0m; I2_0];

x=linspace(0,200,10000)';  

figure(3); 
[tt,yy] = ode45(@(t,y) group_model(t,y,para_m),x,y0_m);
yy12=(yy(:,3)+yy(:,7))*mu1;
yys=cumsum(yy12);

plot(data.ydata(:,1),data.ydata(:,2),'s',tt,yys,'-')
xlabel('Time/day')
ylabel('\mu E(t)')
legend({'Data','\mu E(t)'},'Location','best')
title('Data and fitted model')

modelfun = @(d,th) group_predfun(x,para_m,y0_m);    
nsample = 200;
out = mcmcpred(results,chain,s2chain,x,modelfun,nsample);

figure(4);
mcmcpredplot(out);  
hold on
plot(data.ydata(:,1),data.ydata(:,2),'s'); 
xlabel('Time/day'); ylabel('\mu E(t)');
hold off
title('Predictive envelopes of the model')








