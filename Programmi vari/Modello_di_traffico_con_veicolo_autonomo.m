close all;
N=input('Insert the number of vehicles\n');
T=input('Insert observation time T\n');
L=input('Insert the length of the road (in meters)\n');
k=floor(N/2)+3;
fprintf('\nThe autonomous vehicle is the number %d\n\n',k);
M=10;
r=L/2/pi;
dist=L/N;
x=0:dist:L;
x=x(1:end-1);
x=x(:);
alpha=input('Insert alpha\n');
beta=input('Insert beta\n');
n=input('Insert number of steps\n');
h=T/n;
[xe,ve,xh,vh]=traffic(N,x,L,alpha,beta,n,h,k,M);
anglese=xe/r;
anglesh=xh/r;
circle_xe=r*cos(anglese);
circle_xh=r*cos(anglesh);
circle_ye=r*sin(anglese);
circle_yh=r*sin(anglesh);
theta = linspace(0,2*pi);
x_circle1 = (r-20)*cos(theta);
y_circle1 = (r-20)*sin(theta);
x_circle2 = (r+20)*cos(theta);
y_circle2 = (r+20)*sin(theta);
%%
figure;
time=0:h:T;
for i=1:n
    subplot(2,2,1);
    plot(x_circle1,y_circle1,'black',x_circle2,y_circle2,'black',[circle_xe(1:k-1,i);circle_xe(k+1:end,i)],[circle_ye(1:k-1,i);circle_ye(k+1:end,i)],'.',circle_xe(k,i),circle_ye(k,i),'b.',MarkerSize=9,LineWidth=2);
    axis('square');
    title('Explicit Euler solution');
    posizione = get(gca, 'Position');
    posizione(2) = posizione(2) - 0.13;
    posizione(4) = posizione(4) * 1.5; 
    set(gca, 'Position', posizione);
    subplot(2,2,2);
    plot(x_circle1,y_circle1,'black',x_circle2,y_circle2,'black',[circle_xh(1:k-1,i);circle_xh(k+1:end,i)],[circle_yh(1:k-1,i);circle_yh(k+1:end,i)],'.',circle_xh(k,i),circle_yh(k,i),'b.',MarkerSize=9,LineWidth=2);
    axis('square');
    title('Heun solution');
    posizione = get(gca, 'Position'); 
    posizione(2) = posizione(2) - 0.13;
    posizione(4) = posizione(4) * 1.5;
    set(gca, 'Position', posizione);
    subplot(2,2,3);
    plot(time(1:i),ve(1,1:i),LineWidth=2);
    title('Velocity in time of vehicle 1 (explicit Euler)');
    posizione = get(gca, 'Position');
    posizione(2) = posizione(2) - 0.02;
    posizione(4) = posizione(4) * 0.8; 
    set(gca, 'Position', posizione);
    axis([0 T min(ve(1,:))-0.005 max(ve(1,:))+0.005])
    subplot(2,2,4);
    plot(time(1:i),vh(1,1:i),LineWidth=2);
    title('Velocity in time of vehicle 1 (Heun)');
    posizione = get(gca, 'Position'); 
    posizione(2) = posizione(2) - 0.02;
    posizione(4) = posizione(4) * 0.8; 
    set(gca, 'Position', posizione);
    axis([0 T min(vh(1,:))-0.005 max(vh(1,:))+0.005])
    pause(0.1);
end

%% Cars
k=menu('Do you want to see the dynamics represented with cars?', ...
    'Yes', ...
    'No');
if(k==1)
    posizione_macchinina1=zeros(size(circle_xh(:,1)));
    posizione_macchinina2=zeros(size(circle_yh(:,1)));
    car_img1 = imread('car_image.png'); 
    car_img_resized1 = imresize(car_img1, [12, 15], 'Method', 'bicubic');
    car_img2 = imread('car_image3.png'); 
    car_img_resized2 = imresize(car_img2, [12, 15], 'Method', 'bicubic');
    [img_height, img_width, ~] = size(car_img_resized1);
    figure;
    for i=1:n
        plot(x_circle1,y_circle1,'black',x_circle2,y_circle2,'black',circle_xh(:,i),circle_yh(:,i),'.',LineWidth=2);
        axis('square');
        hold on;
        posizione_macchinina1=[circle_xh(1:k-1,i);circle_xh(k+1:end,i)]-img_width/2;
        posizione_macchinina2=[circle_yh(1:k-1,i);circle_yh(k+1:end,i)]-img_height/2;
        pos1=circle_xh(k,i)-img_width/2;
        pos2=circle_yh(k,i)-img_height/2;
        for j=1:length(posizione_macchinina1)
            if(j~=k)
                image(posizione_macchinina1(j), posizione_macchinina2(j), car_img_resized1);
            else
                image(posizione_macchinina1(j), posizione_macchinina2(j), car_img_resized2);
            end
        end
        hold off;
        title('Traffic dynamic');
        pause(0.1);
    end
end




%%

function [xe,ve,xh,vh]=traffic(N,x,L,alpha,beta,n,h,aut,M)
k=menu('Choose the parameters',...
    'Default: V1=6.75m/s, V2=7.91m/s, C1=0.13m^-1, C2=1.57, l=5m, safety distance=0', ...
    'Customized');
switch k
    case 1
        v1=6.75;
        v2=7.91;
        c1=0.13;
        c2=1.57;
        l=5;
        ds=0;
    case 2
        v1=input('Insert V1 (m/s)\n');
        v2=input('Insert V2 (m/s)\n');
        c1=input('Insert C1 (m^-1)\n');
        c2=input('Insert V2 (constant)\n');
        l=input('Insert l (m)\n');
        ds=input('Insert the safety distance (m)\n');
end
v=L/N*ones(size(x));
v=V(v-l,v1,v2,c1,c2);
gamma=1;
k=menu('Choose the boundary conditions', ...
    'Free flow', ...
    'Periodic');
switch k
    case 1
        v(end+1)=input('Insert the velocity of the ghost vehicle (m/s)\n');
        x(end+1)=x(end)+(x(2)-x(1));
    case 2
        v(end+1)=v(1);
        x(end+1)=x(end)+(x(2)-x(1));
end
g=@f;
pert=menu('Do you want to introduce a perturbation at t=100?', ...
    'Yes', ...
    'No');
[xe,ve]=exp_Euler(h,n,g,x,v,beta,gamma,alpha,v1,v2,c1,c2,l,ds,k,pert,aut,M);
[xh,vh]=Heun(h,n,g,x,v,beta,gamma,alpha,v1,v2,c1,c2,l,ds,k,pert,aut,M);
end






function [v]=f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v,ds,aut,M)          %we are considerating the length of the vehicles and the safety distance
if (aut<length(v))
    app=1/M*sum(dx(aut:aut+M-1)-(l-ds));
    v(1:aut-1)=beta*dv(1:aut-1)./((dx(1:aut-1)-(l-ds)).^(1+gamma))+alpha*(V(dx(1:aut-1)-(l-ds),v1,v2,c1,c2)-v(1:aut-1));
    v(aut)=beta*dv(aut)./((dx(aut)-(l-ds)).^(1+gamma))+alpha*(V(app,v1,v2,c1,c2)-v(aut));
    v(aut+1:end)=beta*dv(aut+1:end)./((dx(aut+1:end)-(l-ds)).^(1+gamma))+alpha*(V(dx(aut+1:end)-(l-ds),v1,v2,c1,c2)-v(aut+1:end));
elseif(aut==length(v))
    app=1/M*sum(dx(aut:aut+M-1)-(l-ds));
    v(1:aut-1)=beta*dv(1:aut-1)./((dx(1:aut-1)-(l-ds)).^(1+gamma))+alpha*(V(dx(1:aut-1)-(l-ds),v1,v2,c1,c2)-v(1:aut-1));
    v(aut)=beta*dv(aut)./((dx(aut)-(l-ds)).^(1+gamma))+alpha*(V(app,v1,v2,c1,c2)-v(aut));
else
    v=beta*dv./((dx-(l-ds)).^(1+gamma))+alpha*(V(dx-(l-ds),v1,v2,c1,c2)-v);
end
end








function [x,v]=Heun(h,n,f,x,v,beta,gamma,alpha,v1,v2,c1,c2,l,ds,k,pert,aut,M) 
dx=x(2:end)-x(1:end-1);
dv=v(2:end)-v(1:end-1);
switch pert
    case 1
        for i=1:99
            v(1:end-1,i+1)=v(1:end-1,i)+h/2*(f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M)+f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M),ds,aut,M));        
            x(1:end-1,i+1)=x(1:end-1,i)+h/2*(v(1:end-1,i)+v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M));
            switch k
                case 1
                    v(end,i+1)=v(end,i);
                    x(end,i+1)=x(end,i)+h*v(end,i);
                case 2
                    v(end,i+1)=v(1,i+1);
                    x(end,i+1)=x(end,i)+h*v(end,i);
            end
            dx=x(2:end,i+1)-x(1:end-1,i+1);
            dv=v(2:end,i+1)-v(1:end-1,i+1);    
        end
        ind=floor(length(x(:,1))/2);
        if (ind>=M)
            v(1:ind-1,101)=v(1:ind-1,100)+h/2*(f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100),ds,aut,M)+f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100)+h*f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100),ds,aut,M),ds,aut,M));    
            v(ind,101)=0.1*v(ind,1);
            v(ind+1:end-1,101)=v(ind+1:end-1,100)+h/2*(f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100),ds,aut,M)+f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100)+h*f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100),ds,aut,M),ds,aut+length(v(ind+1:end-1,101)),M));    %per non avere problemi su M in f
        else
            v(1:ind-1,101)=v(1:ind-1,100)+h/2*(f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100),ds,aut,M)+f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100)+h*f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100),ds,aut,M),ds,aut,M));    
            v(ind,101)=0.1*v(ind,1);
            v(ind+1:end-1,101)=v(ind+1:end-1,100)+h/2*(f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100),ds,aut,M)+f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100)+h*f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100),ds,aut,M),ds,aut-ind,M));
        end
        x(1:end-1,101)=x(1:end-1,100)+h/2*(v(1:end-1,100)+v(1:end-1,100)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,100),ds,aut,M));
        switch k
            case 1
                v(end,101)=v(end,100);
                x(end,101)=x(end,100)+h*v(end,100);
            case 2
                v(end,101)=v(1,101);
                x(end,101)=x(end,100)+h*v(end,100);
        end
        dx=x(2:end,101)-x(1:end-1,101);
        dv=v(2:end,101)-v(1:end-1,101); 
        for i=101:n
            v(1:end-1,i+1)=v(1:end-1,i)+h/2*(f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M)+f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M),ds,aut,M));        
            x(1:end-1,i+1)=x(1:end-1,i)+h/2*(v(1:end-1,i)+v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M));
            switch k
                case 1
                    v(end,i+1)=v(end,i);
                    x(end,i+1)=x(end,i)+h*v(end,i);
                case 2
                    v(end,i+1)=v(1,i+1);
                    x(end,i+1)=x(end,i)+h*v(end,i);
            end
            dx=x(2:end,i+1)-x(1:end-1,i+1);
            dv=v(2:end,i+1)-v(1:end-1,i+1);    
        end
    case 2
        for i=1:n
            v(1:end-1,i+1)=v(1:end-1,i)+h/2*(f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M)+f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M),ds,aut,M));        
            x(1:end-1,i+1)=x(1:end-1,i)+h/2*(v(1:end-1,i)+v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M));
            switch k
                case 1
                    v(end,i+1)=v(end,i);
                    x(end,i+1)=x(end,i)+h*v(end,i);
                case 2                    
                    v(end,i+1)=v(1,i+1);
                    x(end,i+1)=x(end,i)+h*v(end,i);
            end
            dx=x(2:end,i+1)-x(1:end-1,i+1);
            dv=v(2:end,i+1)-v(1:end-1,i+1);    
        end
end
end





function [x,v]=exp_Euler(h,n,f,x,v,beta,gamma,alpha,v1,v2,c1,c2,l,ds,k,pert,aut,M)
dx=x(2:end)-x(1:end-1);
dv=v(2:end)-v(1:end-1);
switch pert
    case 1
        for i=1:99
            v(1:end-1,i+1)=v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M);
            x(1:end-1,i+1)=x(1:end-1,i)+h*v(1:end-1,i);
            switch k
                case 1
                    v(end,i+1)=v(end,i);
                    x(end,i+1)=x(end,i)+h*v(end,i);
                case 2
                    v(end,i+1)=v(1,i+1);
                    x(end,i+1)=x(end,i)+h*v(end,i);
            end
            dx=x(2:end,i+1)-x(1:end-1,i+1);
            dv=v(2:end,i+1)-v(1:end-1,i+1); 
        end        
        ind=floor(length(x(:,1))/2);
        if (ind>=M)
            v(1:ind-1,101)=v(1:ind-1,100)+h*f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100),ds,aut,M);
            v(ind,101)=0.1*v(ind,1);
            v(ind+1:end-1,101)=v(ind+1:end-1,100)+h*f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100),ds,aut+length(v(ind+1:end-1,101)),M);
        else
            v(1:ind-1,101)=v(1:ind-1,100)+h*f(beta,dv(1:ind-1),dx(1:ind-1),gamma,alpha,v1,v2,c1,c2,l,v(1:ind-1,100),ds,aut,M);
            v(ind,101)=0.1*v(ind,1);
            v(ind+1:end-1,101)=v(ind+1:end-1,100)+h*f(beta,dv(ind+1:end),dx(ind+1:end),gamma,alpha,v1,v2,c1,c2,l,v(ind+1:end-1,100),ds,aut-ind,M);
        end
        x(1:end-1,101)=x(1:end-1,100)+h*v(1:end-1,100);
        switch k
            case 1
                v(end,101)=v(end,100);
                x(end,101)=x(end,100)+h*v(end,100);
            case 2
                v(end,101)=v(1,101);
                x(end,101)=x(end,100)+h*v(end,100);
        end
        dx=x(2:end,101)-x(1:end-1,101);
        dv=v(2:end,101)-v(1:end-1,101); 
        for i=101:n
            v(1:end-1,i+1)=v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M);
            x(1:end-1,i+1)=x(1:end-1,i)+h*v(1:end-1,i);
            switch k
                case 1
                    v(end,i+1)=v(end,i);
                    x(end,i+1)=x(end,i)+h*v(end,i);
                case 2
                    v(end,i+1)=v(1,i+1);
                    x(end,i+1)=x(end,i)+h*v(end,i);
            end
            dx=x(2:end,i+1)-x(1:end-1,i+1);
            dv=v(2:end,i+1)-v(1:end-1,i+1); 
        end
    case 2
        for i=1:n
            v(1:end-1,i+1)=v(1:end-1,i)+h*f(beta,dv,dx,gamma,alpha,v1,v2,c1,c2,l,v(1:end-1,i),ds,aut,M);
            x(1:end-1,i+1)=x(1:end-1,i)+h*v(1:end-1,i);
            switch k
                case 1
                    v(end,i+1)=v(end,i);
                    x(end,i+1)=x(end,i)+h*v(end,i);
                case 2
                    v(end,i+1)=v(1,i+1);
                    x(end,i+1)=x(end,i)+h*v(end,i);
            end
            dx=x(2:end,i+1)-x(1:end-1,i+1);
            dv=v(2:end,i+1)-v(1:end-1,i+1); 
        end
end
end



function [dx]=V(dx,v1,v2,c1,c2)           %array form of the function V
dx=v1+v2*tanh(c1*(dx)-c2);
dx=max(0,dx);
end


