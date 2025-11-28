clear;
close all;


N=10000;
R=1;
% n1=[1:4,0.5,1];     
% n2=[1:4,1,0.5];
n1=4;
n2=4;
a=-n1*R;            % Omega=[a,b]
b=n2*R;
p1=1/2;
q1=1/2;
p2=1/2;
q2=1/2;
dv=0.1;
v=a:dv:b;

V=a+(b-a).*rand(N,1);    

T=90;
dt=0.6;
t=0:dt:T;
time=0;
f=1/N*sum(phi(dv,v-V(:,end),2));



while(time+dt<=T)
    V(:,end+1)=V(:,end);
    for i=1:N
        p=rand;
        if(p<dt)
            j=randi(N);
            if(abs(V(i,end-1)-V(j,end-1))<=R)
                p=rand;
                if(p<1/2)
                    V(i,end)=p1*V(i,end-1)+q1*V(j,end-1);
                else
                    V(i,end)=p2*V(i,end-1)+q2*V(j,end-1);
                end
            end
        end
    end
    f(end+1,:)=1/N*sum(phi(dv,v-V(:,end),2));
    time=time+dt;
end
m1=1/N*sum(V);
m2=1/N*sum((V-m1).^2);


%%



% Dimensioni della matrice V
[num_people, num_time_steps] = size(V);

% Definire il range di opinioni con maggiore risoluzione
y = linspace(min(V(:)), max(V(:)), 200); % 200 intervalli per le opinioni

% Definire il tempo
t = linspace(1, num_time_steps, num_time_steps); % 151 intervalli di tempo

% Calcola la densità delle opinioni per ogni intervallo di tempo
density = zeros(length(y), num_time_steps);

for i = 1:num_time_steps
    % Usa una funzione di stima della densità con kernel
    [g, xi] = ksdensity(V(:,i), y, 'Bandwidth', 0.03); % Aggiusta il parametro 'Bandwidth' se necessario
    density(:,i) = g / max(g); % Normalizza per avere un massimo di 1
end


figure;
%for i=1:length(f(:,1))
i=length(f(:,1));
subplot(2,3,3)
    plot(v,f(i,:),LineWidth=1.5);
    title('Distribution of opinions');
    axis([a,b,min(f(:)),max(f(:))+0.2]);
    subplot(2,3,6);
    plot(t,m1,t,m2, LineWidth=1.5);
title('Evolution of the moments')
legend('m_1','m_2',Location='best');
axis([-inf inf -inf inf]);
subplot(2,3,[1,2,4,5]);
contourf(t, y, density, 'LineStyle', 'none');
%colormap('hot');
colorbar;
xlabel('Time');
ylabel('Opinions');
title('Density of opinions over time');
    pause(0.1);
%end



%%

function [p]=phi(dx,x,fun)
switch fun
    case 1
        p=zeros(size(x));
        p(abs(x)<=dx/2)=1/dx;
    case 2
        p=1/dx*B2(x/dx);
end
end


function [b]=B2(x)
b=zeros(size(x));
b(abs(x)<=0.5)=3/4-abs(x(abs(x)<=0.5)).^2;
ind=logical((abs(x)>0.5).*(abs(x)<1.5));
b(ind)=((abs(x(ind))-3/2).^2)/2;
end










