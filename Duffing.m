function x_p=Duffing(t,x,c,k,kappa2,kappa3,eps,sc)
f_nlin=kappa2*x(1)^2+kappa3*x(1)^3-eps*x(2)^3;
x_p=zeros(6,1);



F=exp(sc*t).*([zeros(1) eye(1);-k -c]*x(1:2)-[zeros(1,1);f_nlin]);
DF=exp(sc*t).*([zeros(1) eye(1);-k -c]-[zeros(1,2);2*kappa2*x(1)+3*kappa3*x(1)^2 -3*eps*x(2)^2]);


x_p(1:2)=F;

tmp=(DF)*[x(3) x(4);x(5) x(6)];
x_p(3:6)=reshape(tmp.',4,1);

end