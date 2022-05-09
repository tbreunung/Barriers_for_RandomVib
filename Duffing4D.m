close all
clear all


M=eye(2);
k1=-1;
k2=0.1;
k3=k1;
c1=3.5;
c2=0;%sqrt(3)*0.03;
K= [k1+k2 -k2;-k2 k2+k3];
C=[c1+c2 -c2;-c2 c1+c2];
kappa1=1;
kappa3=1;



 
 
%%
A=[zeros(2) eye(2);-inv(M)*K -inv(M)*C];
[vA,lA]=eig(A);
 
%time parameters
t0=0;
t1=-8;
T_step=-0.01;
tvec=t0:T_step:t1;

%Forcing Direction
B=[0 0 ;0 0; 1 0; 0 1];


% Find zeros
FS=zeros(4,4);
zero1=sqrt(-kappa1/k1);
zero2=sqrt(-kappa3/k3);
%FS(1,:)=zero1.*[1; 0; 0;0];
%FS(2,:)=-FS(1,:);
%FS(3,:)=zero2.*[0; 1; 0;0];
%FS(4,:)=-FS(3,:);
FS(1,:)=[zero1; zero2; 0;0];
FS(2,:)=-FS(1,:);
FS(3,:)=[zero1; -zero2; 0;0];
FS(4,:)=-FS(3,:);
for ii=1:4
 IC=[FS(ii,:)]; %

    [t,tmp]=ode45(@(t,x) LC_example_novar(t,x,M,C,K,kappa1,kappa3),[0 100],IC);
    FSn(ii,:)=tmp(end,1:4);
end

Nphi=100;
phi=linspace(0,2*pi,Nphi);
tol=10^-2;
%p_normal=null([p1 p2].');
 f1=figure;
hold on
f2=figure;
hold on
tra_data=[];
for ii=1:Nphi
    
    x0=tol.*(vA(:,1)*sin(phi(ii))+vA(:,2)*cos(phi(ii)));

    [t,xt]=ode45(@(t,x) LC_example_novar(t,x,M,C,K,kappa1,kappa3),[0 50],x0.');
    
    tra_data=[tra_data; xt(:,1:4) ];   
    figure(f1)
    pl(1,:)=plot3(xt(:,1),xt(:,2),xt(:,3),'m');
    figure(f2)
    pl(1,:)=plot3(xt(:,1),xt(:,2),xt(:,4),'m');
end

[X1, X2]=meshgrid(-1.5:0.03:1.5,-1.5:0.03:1.5);

X3=griddata(tra_data(:,1),tra_data(:,2),tra_data(:,3),X1,X2);
X4=griddata(tra_data(:,1),tra_data(:,2),tra_data(:,4),X1,X2);
figure(f1)
pl(3,:)=surf(X1,X2,X3,'Edgecolor','none','FaceAlpha',0.7,'FaceColor',[0.9290, 0.6940, 0.1250]);

figure(f2)
pl(3,:)=surf(X1,X2,X4,'Edgecolor','none','FaceAlpha',0.7,'FaceColor',[0.9290, 0.6940, 0.1250]);

for ii=1:4
    figure(f1)
    pl(2,:)=plot3(FSn(ii,1),FSn(ii,2),FSn(ii,3),'sk');
    figure(f2)
    pl(2,:)=plot3(FSn(ii,1),FSn(ii,2),FSn(ii,4),'sk');
end
figure(f1)

 xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$\dot{x}_1$','interpreter','latex')

 figure
Datac1=contour( X1,X3,X2,[0 0]);
close
 
%pl(5,:)=plot3(Datac(1,2:end),0.*Datac(1,2:end),Datac(2,2:end),'g','Linewidth',2);
plot3(Datac1(1,2:end),zeros(size(Datac1(1,2:end))),Datac1(2,2:end),'r','Linewidth',2);

 legend(pl,'Trajectories','Stabe Fixed Points','Fitted Slow manifold','location','SouthEast')
axis([-1.2 1.2 -1.5 1.5 -0.2 0.2])

figure(f2)

 xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$\dot{x}_2$','interpreter','latex')

figure
Datac2=contour( X1,X4,X2,[0 0]);
close
 
plot3(Datac2(1,2:end),zeros(size(Datac2(1,2:end))),Datac2(2,2:end),'r','Linewidth',2);

 legend(pl,'Trajectories','Stabe Fixed Points','Fitted Slow manifold','location','SouthEast')
axis([-1.2 1.2 -1.5 1.5 -0.2 0.2])


%%

x0=0.*FS(1,:);



dt=0.1;
Tp=10;
Np=1;
sig=1;
tic
F=@(t,x) Duf_sde(t,x,C,K,kappa1,kappa3);
B=@(t,x) repmat([0.00 ; 0.00 ; 0.00 ;0.01],Np,1);
N=1000000;
parfor ii=1:N
X0=repmat(x0,Np,1)+sig*randn(4*Np,1);
SDE=sde(F,B,'StartState',X0);
[X(:,:,ii),T]=simulate(SDE,Tp/dt,'DeltaTime',dt);
if floor(ii/10000)*10000==ii
disp(['Current progress: ' num2str(round(ii/N *100,2)) '%'])
end

end
T=linspace(0,Tp,Tp/dt);
toc
 %%
Ns=100.*[1 1 1];
bx=1;
x1_min= x0(1)-bx;
x1_max= x0(1)+bx;
x3_min= x0(3)-bx;
x3_max= x0(3)+bx;
x4_min= x0(4)-bx;
x4_max= x0(4)+bx;

tic
x1pts=linspace(x1_min,x1_max,Ns(1));
x3pts=linspace(x3_min,x3_max,Ns(2));
x4pts=linspace(x4_min,x4_max,Ns(3));

 count=zeros(Ns);
% Xs=X;
% idxs=[];
% for xx=1:length(xpts)
%     for yy=1:length(ypts)
%     tmp=Xs-repmat([xpts(xx) ypts(yy)],length(Xs(:,1)),1);
%     idx1=find(abs(tmp(:,1))<dX/2);
%     idx2=find(abs(tmp(idx1,2))<dY/2);
%    Xs(idx1(idx2),:)=[];
%     count(xx,yy)=length(idx2);
%     end
%     xx
% end
% toc
dX1=x1pts(2)-x1pts(1);
dX3=x3pts(2)-x3pts(1);
dX4=x4pts(2)-x4pts(1);

tic
for NN=1:N
    %   tt=length(T)+1;
    for tt=1:length(T)+1
        if abs(X(tt,2,NN)) < dX1/2
            x1_idx=floor((X(tt,1,NN)-x1_min+dX1/2)/dX1)+1;
            
            x3_idx=floor((X(tt,3,NN)-x3_min+dX3/2)/dX3)+1;
            x4_idx=floor((X(tt,4,NN)-x4_min+dX4/2)/dX4)+1;
            
           
            %     xx=floor((X(tt,NN)-x_min+dX/2)/dX)+1;
            %     yy=floor((X(tt,NN+N)-y_min+dX/2)/dY)+1;
            
            if x1_idx<=Ns(1) && x1_idx>0 && x3_idx<=Ns(2)  && x3_idx>0 && x4_idx<=Ns(3) && x4_idx>0 
                count(x1_idx,x3_idx,x4_idx)=count(x1_idx,x3_idx,x4_idx)+1;
            end
        end
    end
end

[vals, idx]=maxk(count(:),200);
[x1s, x3s,x4s]=ind2sub(Ns,idx);

figure 
plot3(Datac1(1,2:end),Datac1(2,2:end),Datac2(2,2:end),'g','LineWidth',4)
%title('Intersection')
xlabel('$q_1$','interpreter','latex','Fontsize',14 )
ylabel('$\dot{q}_1$','interpreter','latex','Fontsize',14 )
zlabel('$\dot{q}_2$','interpreter','latex','Fontsize',14 )

grid on
hold on
%plot3(x1pts(x1s),x3pts(x3s),x4pts(x4s),'sk')
vert0 = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
vert0=repmat([dX1 dX3 dX4],8,1).*(vert0-0.5*ones(size(vert0)));
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
material shiny 

for jj=1:length(x1s)
    cpoint=[x1pts(x1s(jj)) x3pts(x3s(jj)) x4pts(x4s(jj))];
    vert=vert0+repmat(cpoint,8,1); 
     patch('Vertices',vert,'Faces',fac,'FaceColor',[ 0.1 0.1 0.8],'FaceAlpha',0.5,'EdgeColor','none','FaceLighting','gouraud','EdgeLighting','gouraud')
    
end
l = light('Position',[-1 -0.2 0.1],'Style','infinite');

 %set(gcf,'renderer','Painters')


%lighting gouraud
%shp1 = alphaShape(x1pts(x1s).',x3pts(x3s).',x4pts(x4s).');

%pc1 = criticalAlpha(shp1,'all-points');
%shp1.Alpha = pc1+0.01;
%figure;
%hold on; 
%plot(shp1,'EdgeAlpha', 0.0, 'FaceColor','r');
%axis([-1 1 -0.1 0.1 -0.04 0.04])
