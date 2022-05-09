function x_p=Duf_sde(t,x,c,k,kappa2,kappa3)

 Np = numel(x)/2;
 x_p = zeros(2*Np,1);

f_nlin=kappa2.*x(1:Np).^2+kappa3.*x(1:Np).^3;
x_p=[zeros(Np) eye(Np);-diag(repmat(k,1,Np)) -diag(repmat(c,1,Np))]*x-[zeros(Np,1);f_nlin];

% tmp=([zeros(1) eye(1);-k -c]-[zeros(1,2);2*kappa2*x(1)+3*kappa3*x(1)^2 0])*[x(3) x(4);x(5) x(6)];
% x_p(3:6)=reshape(tmp.',4,1);

end