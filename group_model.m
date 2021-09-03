%This is the covid model with two groups

function ydot = group_model(t,y,k)

global mu1 mu2

S1=y(1); V=y(2); E1=y(3); U1=y(4); I1=y(5); S2=y(6); E2=y(7); U2=y(8); I2=y(9); 

beta11=k(1); beta12=k(2); omega=k(3); 
beta11v=k(4); epsilon=k(5); beta12v=k(6); 
mu1=1/3; rho=0.6; 
eta11=1/7; eta12=eta11; tau1=1/4.6; delta1=k(7); beta21=k(8);
beta22=k(9); mu2=mu1; eta21=eta11; eta22=eta21;
tau2=tau1; delta2=k(10); 

Delta11=beta11*U1*S1;
Delta12=beta12*U2*S1;
Delta11v=beta11v*(1-epsilon)*U1*V;
Delta12v=beta12v*(1-epsilon)*U2*V;
Delta21=beta21*U2*S2;
Delta22=beta22*U1*S2;
 

   ydot = [-Delta11-Delta12-omega*S1;
       omega*S1-Delta11v-Delta12v;
       Delta11+Delta12+Delta11v+Delta12v-mu1*E1;
       mu1*rho*E1-eta11*U1-tau1*U1;
       mu1*(1-rho)*E1+tau1*U1-eta12*I1-delta1*I1;
       -Delta21-Delta22;
       Delta21+Delta22-mu2*E2;
       mu2*rho*E2-eta21*U2-tau2*U2;
       mu2*(1-rho)*E2+tau2*U2-eta22*I2-delta2*I2];




