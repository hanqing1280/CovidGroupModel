% This function returns the solution E1(t), E2(t) 

function  ys= group_predfun(time,k,y0)

[~,yp] = ode45(@(t,y) group_model(t,y,k),time,y0);
m=length(time);
ys1=zeros(m,1);
ys2=zeros(m,1);
ys1(:)=yp(:,3); % E_1(t)
ys2(:)=yp(:,7); % E_2(t)
ys12=ys1+ys2;
ys=cumsum(ys12);




