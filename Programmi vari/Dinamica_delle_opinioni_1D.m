close all;
clear;

%transient= time delay for reaching the equilibrium

rng default;
N=80;
n=200;
T=20;
dt=T/n;
epsilon=0.25;            %USE THIS PARAMETER TO CHANGE INTERACTION RANGE
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
mean=1/N*sum(x);
variance=sum(x.^2)/N-mean.^2;
if(length(x(:,1))==N)
    N1=sum(abs(x-x(1,:))<epsilon);
elseif(length(x(:,1))==N+1)
    N1=sum(abs(x(1:end-1,:)-x(1,:))<epsilon);
end

%for the transient
toll=0.0000001;
eq_time=find(abs(variance(2:end)-variance(1:end-1))<toll,1);      
eq_time=eq_time*dt;

algorithm_duration=toc;     %tic is inside the Hegselmann_Krause function; how much does it change when the number of agents changes?

figure;
subplot(2,2,[1,2]);
plot(0:dt:T,x);
title('Opinion dynamics');
subplot(2,2,3);
plot(0:dt:T,mean,0:dt:T,variance);
title('Mean and variance of opinions');
legend('Mean','Variance',Location='best');
subplot(2,2,4);
plot(0:dt:T,N1);
title('Variation in time of N_1');
if strcmp(type,'Non-symmetric')
    subplot(2,2,4);
    plot(0:dt:T,multiplicity);
    title('Multiplicity of \lambda=1');
end



        


%%


function [out]=Hegselmann_Krause(x,T,dt,epsilon,type)
k=menu('Choose the type of model', ...
        'Without leader', ...
        'With leader');
switch k
    case 1
        tic;        %for the study duration
        fun=@matrix_f;
        n=T/dt;
        out=Heun(0,dt,n,fun,x,epsilon,type);
    case 2
        tic;        %for the study duration
        fun=@matrix_f_leader;
        x(end+1)=0.7;
        n=T/dt;
        out=Heun_leader(0,dt,n,fun,x,epsilon,type);
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


function [x]=initial_value(N)
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
end



function [out]=matrix_f(~,u,epsilon,type)
N=length(u);
switch type
    case 'Symmetric'
        K=zeros(N);
        interacting=abs(u-u')<epsilon;        %interacting(i,j) is 1 if i and j interact, 0 otherwise
        Ni=sum(interacting);
        K(interacting)=1/N;
        K=K-diag(Ni/N);                     %if epsilon>max||xj-xi|| then Ni=N and the interaction is global
        u=K*u;
        out=u;                                      %I use out because I want a different amount of outputs depending on cases
        
    case 'Non-symmetric'
        K=zeros(N);
        interacting=abs(u-u')<epsilon;
        Ni=sum(interacting);
        Ni=Ni(:);
        Ni=repmat(Ni,1,N);
        K(interacting)=1./Ni(interacting);                  %this is the matrix A from theory
        multiplicity=sum(abs(eig(K)-1)<0.00000000000001);
        K=K-eye(N);
        u=K*u;
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




function [out]=matrix_f_leader(~,u,epsilon,type)                %u is N+1 long
N=length(u)-1;
switch type
    case 'Symmetric'
        K=zeros(N);
        interacting=abs(u(1:N)-u(1:N)')<epsilon;        %interacting(i,j) is 1 if i and j interact, 0 otherwise
        Ni=sum(interacting);
        K(interacting)=1/N;
        K=K-diag(Ni/N);                     %if epsilon>max||xj-xi|| then Ni=N and the interaction is global
        u(1:N)=K*u(1:N)+(u(end)-u(1:N));
        out=u(1:N);                                      %u(end) is constant
        
    case 'Non-symmetric'
        K=zeros(N);
        interacting=abs(u(1:N)-u(1:N)')<epsilon;
        Ni=sum(interacting);
        Ni=Ni(:);
        Ni=repmat(Ni,1,N);
        K(interacting)=1./Ni(interacting);                  %this is the matrix A from theory
        multiplicity=sum(abs(eig(K)-1)<0.00000000000001);
        K=K-eye(N);
        u(1:N)=K*u(1:N)+(u(end)-u(1:N));
        out={u(1:N),multiplicity};
end
end



function [out]=Heun_leader(t0,h,n,f,u,e,t)                     %e and t are needed for the function f; we split in two cases in order to save the multiplicity
switch t
    case 'Symmetric'
        for i=1:n
            u(1:end-1,i+1)=u(1:end-1,i)+h/2*(f(t0+(i-1)*h,u(:,i),e,t)+f(t0+i*h,u(:,i)+h*[f(t0+(i-1)*h,u(:,i),e,t);0],e,t));      %I use the array for the count otherwise the leader would be ignored
            u(end,i+1)=u(end,i);
        end
        out=u;
    case 'Non-symmetric'
        multiplicity=zeros(1,n+1);
        for i=1:n
            result=f(t0+(i-1)*h,u(:,i),e,t);                %the output is a cell array
            first=result{1};
            multiplicity(i)=result{2};
            result=f(t0+i*h,u(:,i)+h*[first;0],e,t);
            second=result{1};
            u(1:end-1,i+1)=u(1:end-1,i)+h/2*(first+second);
            u(end,i+1)=u(end,i);
        end
        result=f(t0+n*h,u(:,n+1),e,t);
        multiplicity(n+1)=result{2};
        out={u,multiplicity};
end
end




