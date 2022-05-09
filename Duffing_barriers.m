close all
clear all


k=-1;
c=3.5;
kappa2=0;
kappa3=1;
%eps=0.01;

x0=sqrt(-k/kappa3);
y0=0;


dt=0.1;
Tp=10;
Np=1;
sig=1;


tic
F=@(t,x) Duf_sde(t,x,c,k,kappa2,kappa3);
B=@(t,x) repmat([0.00; 0.1],Np,1);

N=10^6;

%N=1000000;
parfor ii=1:N
    X0=repmat([x0;y0],Np,1)+sig*randn(2*Np,1);
    SDE=sde(F,B,'StartState',X0);
    [X(:,:,ii),T]=simulate(SDE,Tp/dt,'DeltaTime',dt);
    if floor(ii/10000)*10000==ii
        disp(['Current progress: ' num2str(round(ii/N *100,2)) '%'])
    end
end
T=linspace(0,Tp,Tp/dt);
toc
% %plot(squeeze(X(:,1,:)),squeeze(X(:,2,:)))
% axis equal
% hold on
% A=[0 1;    2*k -c];
% [vA,lA]=eig(A);
% plot([x0 x0]+0.25.*vA(1,1).*[-1 1],[y0 y0]+0.25.*vA(2,1).*[-1 1],'-g');
% plot([x0 x0]+0.25.*vA(1,2).*[-1 1],[y0 y0]+0.25.*vA(2,2).*[-1 1],'-r');
% 
%%
Nx=100;
Ny=100;
bx=0.5;
x_min=x0-bx;
x_max=x0+bx;
y_min=y0-bx;
y_max=y0+bx;
tic
xpts=linspace(x_min,x_max,Nx);
ypts=linspace(y_min,y_max,Ny);
count=zeros(length(xpts),length(ypts));
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
dX=xpts(2)-xpts(1);
dY=ypts(2)-ypts(1);
tic
for NN=1:N
    %   tt=length(T)+1;
    for tt=1:round(length(T)/10)+1
        xx=floor((X(tt,1,NN)-x_min+dX/2)/dX)+1;
        yy=floor((X(tt,2,NN)-y_min+dX/2)/dY)+1;
        %     xx=floor((X(tt,NN)-x_min+dX/2)/dX)+1;
        %     yy=floor((X(tt,NN+N)-y_min+dX/2)/dY)+1;
        
        if xx<=Ny && xx>0 && yy<=Nx  && yy>0
            count(xx,yy)=count(xx,yy)+1;
        end
        
    end
end
toc


A=[0 1;    2*k -c];
[vA,lA]=eig(A);

 [SSM1, SSM2]=Duf_SSM(lA(1,1),lA(2,2),k,kappa2,kappa3,1);
 SSM1=SSM1+repmat([x0;y0],1,100);
 SSM2=SSM2+repmat([x0;y0],1,100);
 
 
figure
imagesc( [x_min x_max],[y_max y_min],flipud( (count)'))

h=colorbar;
hold on
% plot([x0 x0]+bx.*vA(1,1).*[-1 1],[y0 y0]+bx.*vA(2,1).*[-1 1],'-g');
% plot([x0 x0]+bx.*vA(1,2).*[-1 1],[y0 y0]+bx.*vA(2,2).*[-1 1],'-r');

%pl=plot(SSM2(1,:), SSM2(2,:) ,'Linewidth',2);
hold on
[~,sep]=ode45(@(t,x) Duffing(t,x,c,k,kappa2,kappa3,0,0),[0 100],[0.01 0 1 0 0 1]);
pl=plot(sep(:,1),sep(:,2),'k');

legend(pl, 'Separatrix','AutoUpdate','off')

ylabel(h, 'count')
set(gca,'YDir','normal') 
xlabel('Postion')
ylabel('Velocity')

figure
 
% cnt_tmp(cnt_tmp>10)=NaN;
imagesc( [x_min x_max],[y_max y_min],flipud((count)'))

h=colorbar;
hold on
% plot([x0 x0]+bx.*vA(1,1).*[-1 1],[y0 y0]+bx.*vA(2,1).*[-1 1],'-g');
% plot([x0 x0]+bx.*vA(1,2).*[-1 1],[y0 y0]+bx.*vA(2,2).*[-1 1],'-r');

pl=plot(SSM2(1,:),SSM2(2,:),'Linewidth',2);
% [~,sep]=ode45(@(t,x) Duffing(t,x,c,k,kappa2,kappa3,0,0),[0 100],[0.01 0 1 0 0 1]);
% pl=plot(sep(:,1),sep(:,2),'g');
pl=plot(SSM2(1,:),SSM2(2,:),'Linewidth',2);
%legend(pl,['SSM tngt. to E with Re(\lambda)=' num2str(real(lA(1,1)))],['SSM tngt. to E with Re(\lambda)=' num2str(real(lA(2,2)))],'Separatrix','AutoUpdate','off')

ylabel(h, 'count')
set(gca,'YDir','normal') 
xlabel('Postion')
ylabel('Velocity')



%%
% figure
% [M,cs]=contour(xpts,ypts,log(count') );
% hold on
% plot([x0 x0]+bx.*vA(1,1).*[-1 1],[y0 y0]+bx.*vA(2,1).*[-1 1],'-g');
% plot([x0 x0]+bx.*vA(1,2).*[-1 1],[y0 y0]+bx.*vA(2,2).*[-1 1],'-r');
%pl=plot(SSM2(1,:),SSM2(2,:),'Linewidth',2);
%plot(sep(:,1),sep(:,2),'k')
axis([x_min x_max y_min y_max])
phi=linspace(pi,2*pi);

%phi=linspace(pi/2,3/2*pi);
R=5;
X0=repmat([1;0],1,R);
%dr=0.1;
dr=0.1;
for rr=1:R
     cnt=zeros(100,1);
for phi_it=1:100
    phi_it
 
 
vec1=diff(X,1);
%vec2=dr.*[cos(phi(phi_it));-sin(phi(phi_it))];
vec2=dr.*[cos(phi(phi_it));-sin(phi(phi_it))];
for nn=1:10:N
  
for tt=1:length(T)-1  %length(T)+1  
    
%[~,Xt]=ode45(@(t,x) Duffing(t,x,c,k,kappa2,kappa3),[0 -dt*tt],[squeeze(X(tt,:,nn))  1 0 0 1]);      
%vec1=squeeze(X(1,:,nn)).'-[Xt(1);Xt(2)];
   M=[vec1(tt,:,nn).' -vec2];
   al=M\(-squeeze(X(tt,:,nn)).'+X0(:,rr));
   
    if al(1)<=1 && al(1)>=0 && al(2)<=1 && al(2)>=0
        cnt(phi_it)=cnt(phi_it)+1;
        
    end
end
end

end


[~,idx]=min(cnt);
X0(:,rr+1)=X0(:,rr)+dr.*[cos(phi(idx));-sin(phi(idx))];
end
pl(2,:)=plot(X0(1,:),X0(2,:),'k','Linewidth',1);
legend(pl,'SSM' ,'Barrier','AutoUpdate','off')
