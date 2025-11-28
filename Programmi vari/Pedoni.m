close all;
clear;

N=100;
m=1;
r=0.025;
yw1=2.1;
yw2=-2.1;
u(1:50,1:2)=[-7+4*rand(N/2,1),-1+2*rand(N/2,1)];
u(1:50,3:4)=[ones(50,1),zeros(50,1)];
u(51:100,1:2)=[3+4*rand(N/2,1),-1+2*rand(N/2,1)];
u(51:100,3:4)=[-ones(50,1),zeros(50,1)];

v0=normrnd(1.34,0.26,N,1);
e0(1:50,1:2)=[ones(50,1),zeros(50,1)];
e0(51:100,1:2)=[-ones(50,1),zeros(50,1)];
tau=0.1;
Ar=2*10^2;
Aw=Ar;
Br=0.08;
Bw=Br;
k=0;
kappa=k;
% T=input('Observation time\n');
% n=input('Number of nodes\n');
T=100;
n=1000;
dt=T/n;
u=heun(u,dt,n,m,v0,e0,tau,20,Br,Aw,Bw,r,k,kappa,yw1,yw2);   


%%
figure;
i=1;
while(any((-10.2<=u(:,1,i)) & (u(:,1,i)<=10.2)) && i<=n)
    plot([-10.2,10.2],[yw1+0.2,yw1+0.2],'black-',[-10.2,10.2],[yw2-0.2,yw2-0.2],'black-',u(1:50,1,i),u(1:50,2,i),'bo',u(51:end,1,i),u(51:end,2,i),'ro',MarkerSize=3,LineWidth=4);
    axis([-10 10 yw2-1 yw1+1]);
    title('Pedestrian dynamics');
    pause(0.1);
    i=i+1;
end




%%




function[fun]=f(u,m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2)
fun(:,1:2)=u(:,3:4);
psyc=fr(u,Ar,r,Br,k,kappa);
sum1=sum(psyc(:,:,1)-diag(diag(psyc(:,:,1))),2);
sum2=sum(psyc(:,:,2)-diag(diag(psyc(:,:,2))),2);
sum3=[sum1,sum2];
fun(:,3:4)=1./m.*fd(u,m,v0,e0,tau)+sum3+fw(u,Aw,r,Bw,yw1)+fw(u,Aw,r,Bw,yw2);
end





function [f]=fd(u,m,v0,e0,tau)
f=m.*(v0.*e0-u(:,3:4))./tau;
end



function [f]=fr(u,A,r,B,k,kappa)   
d=pdist(u(:,1:2));
d=squareform(d);
diff=(r+r')-d;
g=max(0,diff);
n=(u(:,1)-u(:,1)')./(d+eye(size(d)));
n(:,:,2)=(u(:,2)-u(:,2)')./(d+eye(size(d)));
t=zeros(size(n));
t(:,:,1)=-n(:,:,2);
t(:,:,2)=n(:,:,1);
dv(:,:,1)=(u(:,3)-u(:,3)')';
dv(:,:,2)=(u(:,4)-u(:,4)')';
dv=sum(dv.*t,3);
f=(A.*exp(diff./B)+k*g).*n+kappa.*g.*dv.*t;
end



function [f]=fw(u,A,r,B,yw)
dw=abs(u(:,2)-yw);
nw=sign(u(:,2)-yw);
nw=[zeros(size(nw)),nw];
f=(A.*exp((r-dw)./B)).*nw;
end



function [u]=heun(u,h,n,m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2)
for i=1:n
    u(:,:,i+1)=u(:,:,i)+h/2*(f(u(:,:,i),m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2)+f(u(:,:,i)+h*f(u(:,:,i),m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2),m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2));
end
end









