function x_p=LC_example_novar(t,x,M,C,K,kappa1,kappa3 )
dim=min(size(M));
f_nlin=[kappa1*x(1)^3 ;kappa3*x(2)^3 ];
x_p=zeros(2*dim ,1);
x_p(1:2*dim)=[zeros(dim) eye(dim);-inv(M)*K -inv(M)*C]*x(1:2*dim)-[zeros(dim,1);f_nlin];


end