close all;
clear;

N=40;
n=70;
T=10;
dt=T/n;
epsilon=0.3;            %USE THIS PARAMETER TO CHANGE INTERACTION RANGE
type=type_of_interactions;
x=initial_value(N);
output=Hegselmann_Krause(x,T,dt,epsilon,type);
switch type
    case 'Symmetric'
        x=output;              
    case 'Non-symmetric'
        x=output{1};
        multiplicity=output{2};
end
if(length(x(:,1))==2*N)
    mean1=1/N*sum(x(1:N,:));                    %we don't have any leader
    variance1=sum(x(1:N,:).^2)/N-mean1.^2;
    mean2=1/N*sum(x(N+1:end,:));
    variance2=sum(x(N+1:end,:).^2)/N-mean2.^2;
    N1=sum((abs(x(1:N,:)-x(1,:))<epsilon) & (abs(x(N+1:end,:)-x(N+1,:))<epsilon));

    if strcmp(type,'Symmetric')
        figure;
        subplot(2,2,1);
        plot(0:dt:T,x(1:N,:));
        title('Opinion 1 dynamics');
        subplot(2,2,2);
        plot(0:dt:T,x(N+1:end,:));
        title('Opinion 2 dynamics');
        subplot(2,2,3);
        plot(0:dt:T,mean1,0:dt:T,variance1,0:dt:T,mean2,0:dt:T,variance2);
        title('Mean and variance of opinions');
        legend('Mean 1','Variance 1','Mean 2','Variance 2',Location='best');
        subplot(2,2,4);
        plot(0:dt:T,N1);
        title('Variation in time of N_1');
    elseif strcmp(type,'Non-symmetric')
        figure;
        subplot(3,2,1);
        plot(0:dt:T,x(1:N,:));
        title('Opinion 1 dynamics');
        subplot(3,2,2);
        plot(0:dt:T,x(N+1:end,:));
        title('Opinion 2 dynamics');
        subplot(3,2,3);
        plot(0:dt:T,mean1,0:dt:T,variance1,0:dt:T,mean2,0:dt:T,variance2);
        title('Mean and variance of opinions');
        legend('Mean 1','Variance 1','Mean 2','Variance 2',Location='best');
        subplot(3,2,4);
        plot(0:dt:T,N1);
        title('Variation in time of N_1');
        subplot(3,2,5);
        plot(0:dt:T,multiplicity);
        title('Multiplicity of \lambda=1');
     end

elseif(length(x(:,1))==2*N+2)
    mean1=1/N*sum([x(1:N,:);x(end-1,:)]);                    %we consider the leader
    variance1=sum([x(1:N,:);x(end-1,:)].^2)/N-mean1.^2;
    mean2=1/N*sum([x(N+1:end-2,:);x(end,:)]);
    variance2=sum([x(N+1:end-2,:);x(end,:)].^2)/N-mean2.^2;
    N1=sum((abs(x(1:N,:)-x(1,:))<epsilon) & (abs(x(N+1:end-2,:)-x(N+1,:))<epsilon));

     if strcmp(type,'Symmetric')
        figure;
        subplot(2,2,1);
        plot(0:dt:T,[x(1:N,:);x(end-1,:)]);
        title('Opinion 1 dynamics');
        subplot(2,2,2);
        plot(0:dt:T,[x(N+1:end-2,:);x(end,:)]);
        title('Opinion 2 dynamics');
        subplot(2,2,3);
        plot(0:dt:T,mean1,0:dt:T,variance1,0:dt:T,mean2,0:dt:T,variance2);
        title('Mean and variance of opinions');
        legend('Mean 1','Variance 1','Mean 2','Variance 2',Location='best');
        subplot(2,2,4);
        plot(0:dt:T,N1);
        title('Variation in time of N_1');
     elseif strcmp(type,'Non-symmetric')
        figure;
        subplot(3,2,1);
        plot(0:dt:T,[x(1:N,:);x(end-1,:)]);
        title('Opinion 1 dynamics');
        subplot(3,2,2);
        plot(0:dt:T,[x(N+1:end-2,:);x(end,:)]);
        title('Opinion 2 dynamics');
        subplot(3,2,3);
        plot(0:dt:T,mean1,0:dt:T,variance1,0:dt:T,mean2,0:dt:T,variance2);
        title('Mean and variance of opinions');
        legend('Mean 1','Variance 1','Mean 2','Variance 2',Location='best');
        subplot(3,2,4);
        plot(0:dt:T,N1);
        title('Variation in time of N_1');
        subplot(3,2,5);
        plot(0:dt:T,multiplicity);
        title('Multiplicity of \lambda=1');
    end
end




        


%%


function [out]=Hegselmann_Krause(x,T,dt,epsilon,type)
k=menu('Choose the type of model', ...
        'Without leader', ...
        'With leader');
switch k
    case 1
        fun=@matrix_f2;
        n=T/dt;
        out=Heun(0,dt,n,fun,x,epsilon,type);
    case 2
        fun=@matrix_f2_leader;
        x(end+1)=0.3;
        x(end+1)=0.4;
        n=T/dt;
        out=Heun2_leader(0,dt,n,fun,x,epsilon,type);
end
end




function [type]=type_of_interactions ()
k=menu('Choose the type of interactions', ...
        'Symmetric', ...
        'Non-symmetric');
    
    switch k
        case 1
            type=char('Symmetric');
        case 2
            type=char('Non-symmetric');
    end
end


function [u]=initial_value(N)
k=menu('Choose x_0 distribution', ...
        'x_0~U(0,1)', ...
        'x_0~N(0.5,0.3^4)', ...
        'x_0~P(0.3,0.15,1,3)');
    
    switch k
        case 1
            x=rand(N,1);
        case 2
            x=normrnd(0.5,0.3^4,N,1);
        case 3
            x=pearsrnd(0.3,0.15,1,3,N,1);
    end
k=menu('Choose y_0 distribution', ...
        'y_0~U(0,1)', ...
        'y_0~N(0.5,0.3^4)', ...
        'y_0~P(0.3,0.15,1,3)');
    
    switch k
        case 1
            y=rand(N,1);
        case 2
            y=normrnd(0.5,0.3^4,N,1);
        case 3
            y=pearsrnd(0.3,0.15,1,3,N,1);
    end
    u=[x;y];
end



function [out]=matrix_f2(~,u,epsilon,type)           %without leader
N=length(u);
switch type
    case 'Symmetric'
        K=zeros(N/2);
        interacting=((abs(u(1:N/2)-u(1:N/2)')<epsilon) & (abs(u(N/2+1:end)-u(N/2+1:end)')<epsilon));      
        Ni=sum(interacting);
        K(interacting)=1/(N/2);
        K=K-diag(Ni/(N/2));               
        u(1:N/2)=K*u(1:N/2);
        u(N/2+1:end)=K*u(N/2+1:end);
        out=u;                                      %I use out because I want a different amount of outputs depending on cases
        
    case 'Non-symmetric'
        K=zeros(N/2);
        interacting=((abs(u(1:N/2)-u(1:N/2)')<epsilon) & (abs(u(N/2+1:end)-u(N/2+1:end)')<epsilon));
        Ni=sum(interacting);
        Ni=Ni(:);
        Ni=repmat(Ni,1,N/2);
        K(interacting)=1./Ni(interacting);                 
        multiplicity=sum(abs(eig(K)-1)<0.00000000000001);
        K=K-eye(N/2);
        u(1:N/2)=K*u(1:N/2);
        u(N/2+1:end)=K*u(N/2+1:end);
        out={u,multiplicity};
end
end




function [out]=Heun(t0,h,n,f,u,e,t)                     %e and t are needed for the function f; we split in two cases in order to save the multiplicity
switch t
    case 'Symmetric'
        for i=1:n
            u(:,i+1)=u(:,i)+h/2*(f(t0+(i-1)*h,u(:,i),e,t)+f(t0+i*h,u(:,i)+h*f(t0+(i-1)*h,u(:,i),e,t),e,t));
        end
        out=u;
    case 'Non-symmetric'
        multiplicity=zeros(1,n+1);
        for i=1:n
            result=f(t0+(i-1)*h,u(:,i),e,t);                %the output is a cell array
            first=result{1};
            multiplicity(i)=result{2};
            result=f(t0+i*h,u(:,i)+h*first,e,t);
            second=result{1};
            u(:,i+1)=u(:,i)+h/2*(first+second);
        end
        result=f(t0+n*h,u(:,n+1),e,t);
        multiplicity(n+1)=result{2};
        out={u,multiplicity};
end
end




function [out]=matrix_f2_leader(~,u,epsilon,type)                %u is N+1 long
N=length(u)-2;
switch type
    case 'Symmetric'
        K=zeros(N/2);
        interacting=((abs(u(1:N/2)-u(1:N/2)')<epsilon) & (abs(u(N/2+1:end-2)-u(N/2+1:end-2)')<epsilon));       
        Ni=sum(interacting);
        K(interacting)=1/(N/2);
        K=K-diag(Ni/(N/2));                    
        u(1:N/2)=K*u(1:N/2)+(u(end-1)-u(1:N/2));
        u(N/2+1:end-2)=K*u(N/2+1:end-2)+(u(end)-u(N/2+1:end-2));
        out=u(1:N);                                      %u(end) is constant
        
    case 'Non-symmetric'
        K=zeros(N/2);
        interacting=((abs(u(1:N/2)-u(1:N/2)')<epsilon) & (abs(u(N/2+1:end-2)-u(N/2+1:end-2)')<epsilon)); 
        Ni=sum(interacting);
        Ni=Ni(:);
        Ni=repmat(Ni,1,N/2);
        K(interacting)=1./Ni(interacting);               
        multiplicity=sum(abs(eig(K)-1)<0.00000000000001);
        K=K-eye(N/2);
        u(1:N/2)=K*u(1:N/2)+(u(end-1)-u(1:N/2));
        u(N/2+1:end-2)=K*u(N/2+1:end-2)+(u(end)-u(N/2+1:end-2));
        out={u(1:N),multiplicity};
end
end



function [out]=Heun2_leader(t0,h,n,f,u,e,t)                     %e and t are needed for the function f; we split in two cases in order to save the multiplicity
switch t
    case 'Symmetric'
        for i=1:n
            u(1:end-2,i+1)=u(1:end-2,i)+h/2*(f(t0+(i-1)*h,u(:,i),e,t)+f(t0+i*h,u(:,i)+h*[f(t0+(i-1)*h,u(:,i),e,t);0;0],e,t));      %I use the array for the count otherwise the leader would be ignored
            u(end-1,i+1)=u(end-1,i);
            u(end,i+1)=u(end,i);
        end
        out=u;
    case 'Non-symmetric'
        multiplicity=zeros(1,n+1);
        for i=1:n
            result=f(t0+(i-1)*h,u(:,i),e,t);                %the output is a cell array
            first=result{1};
            multiplicity(i)=result{2};
            result=f(t0+i*h,u(:,i)+h*[first;0;0],e,t);
            second=result{1};
            u(1:end-2,i+1)=u(1:end-2,i)+h/2*(first+second);
            u(end-1,i+1)=u(end-1,i);
            u(end,i+1)=u(end,i);
        end
        result=f(t0+n*h,u(:,n+1),e,t);
        multiplicity(n+1)=result{2};
        out={u,multiplicity};
end
end