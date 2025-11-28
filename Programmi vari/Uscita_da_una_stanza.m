close all;
clear;
N=100;
m=1;
r=0.025;
yw1=2.1;
yw2=-2.1;
xw1=-10;
xw2=10;
xdoor=0;
ydoor=yw1;
doorwidth=r+0.2;
u(:,1:2)=[-5+10*rand(N,1),-1+2*rand(N,1)];
u(:,3:4)=[zeros(N,1),ones(N,1)];

v0=normrnd(1.34,0.26,N,1);
e0(:,1:2)=[-u(:,1),yw1*5-u(:,2)]./(u(:,1).^2+(yw1*5-u(:,2)).^2).^(1/2);
tau=0.05;
Ar=10;
Aw=2*10^2;
Br=0.08;
Bw=Br;
k=1.2;
kappa=2.4;
T=100;
n=1000;
dt=T/n;
u=heun(u,dt,n,m,v0,e0,tau,10,Br,Aw,Bw,r,k,kappa,yw1,yw2,xw1,xw2,doorwidth);


%%
figure;
for i=1:n
    plot(u(:,1,i),u(:,2,i),'bo',MarkerSize=3,LineWidth=4); 
    axis([xw1-1 xw2+1 yw2-1 yw1+1]);
    hold on;
    rectangle('Position',[xw1-0.1 yw2-0.1 xw2-xw1+0.2 yw1-yw2+0.2], LineWidth=4);
    plot([-doorwidth-0.1,doorwidth+0.1],[yw1+0.1,yw1+0.1],'w-',LineWidth=4);
    hold off;
    title('Pedestrian dynamics');
    pause(0.1);
end




%%




function[fun]=f(u,m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2,xw1,xw2,doorwidth)
fun(:,1:2)=u(:,3:4);
psyc=fr(u,Ar,r,Br,k,kappa);
sum1=sum(psyc(:,:,1)-diag(diag(psyc(:,:,1))),2);
sum2=sum(psyc(:,:,2)-diag(diag(psyc(:,:,2))),2);
sum3=[sum1,sum2];
fun(:,3:4)=1./m.*fd(u,m,v0,e0,tau)+sum3+fwdoor(u,Aw,r,Bw,yw1,doorwidth)+fwy(u,Aw,r,Bw,yw2)+fwx(u,Aw,r,Bw,xw1)+fwx(u,Aw,r,Bw,xw2);
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



function [f]=fwy(u,A,r,B,yw)
dw=abs(u(:,2)-yw);
nw=-sign(yw-u(:,2));
nw=[zeros(size(nw)),nw];
f=(A.*exp((r-dw)./B)).*nw;
end



function [f]=fwdoor(u,A,r,B,yw,doorwidth)
dw=abs(u(:,2)-yw);
nw=-sign(yw-u(:,2));
nw=[zeros(size(nw)),nw];
f=(A.*exp((r-dw)./B)).*nw;
ind=abs(u(:,1))<doorwidth;
f(ind,:)=0;
end


function [f]=fwx(u,A,r,B,xw)
dw=abs(u(:,1)-xw);
nw=-sign(xw-u(:,1));
nw=[nw,zeros(size(nw))];
f=(A.*exp((r-dw)./B)).*nw;
end



function [u]=heun(u,h,n,m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2,xw1,xw2,doorwidth)
for i=1:n
    u(:,:,i+1)=u(:,:,i)+h/2*(f(u(:,:,i),m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2,xw1,xw2,doorwidth)+f(u(:,:,i)+h*f(u(:,:,i),m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2,xw1,xw2,doorwidth),m,v0,e0,tau,Ar,Br,Aw,Bw,r,k,kappa,yw1,yw2,xw1,xw2,doorwidth));
    e0(:,1:2)=[-u(:,1),yw1*3-u(:,2)]./(u(:,1).^2+(yw1*5-u(:,2)).^2).^(1/2);   
end
end

