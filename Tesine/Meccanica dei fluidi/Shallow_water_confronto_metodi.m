clc;
close all; clear;


% togliere il commento ad una delle sezioni di fila per le simulazioni preimpostate
% è presente un waitforbuttonpress
% si aprono molteplici schede per i grafici


%LAR
B=@(x) 0.2+0.5.*exp(1-1./(1-((abs(x-10)./5).^2))).*(x<15).*(x>5)+0*x;
% B=@(x) 0.4+0.1.*sin(pi.*x)+0.002.*x.^2;
h0=@(x) 2.2-B(x)+0.05.*(x<4.5).*(x>3.8);
regime='lar';
% 
% % super (cond a sx)
% B=@(x) 0.2+0.5.*exp(1-1./(1-((abs(x-10)./5).^2))).*(x<15).*(x>5)+0*x;
% h0=@(x) 2.2-B(x);
% qcond=24;
% hcond=2;
% regime='super';
% 
% % sub (cond q a sx e h a dx)
% B=@(x) 0.2+0.5.*exp(1-1./(1-((abs(x-10)./5).^2))).*(x<15).*(x>5)+0*x;
% h0=@(x) 2.2-B(x);
% qcond=3.42;
% hcond=2;
% regime='sub';
% 
% %trans1 (cond q a sx e h a dx sopra una soglia)
% B=@(x) 0.2+0.5.*exp(1-1./(1-((abs(x-10)./5).^2))).*(x<15).*(x>5)+0*x;
% h0=@(x) 0.8-B(x);
% qcond=1.58;
% hcond=0.6;
% regime='trans';
% 
% %trans2 (cond q a sx e h a dx sopra una soglia)
% B=@(x) 0.2+0.5.*exp(1-1./(1-((abs(x-10)./5).^2))).*(x<15).*(x>5)+0*x;
% h0=@(x) 0.8-B(x);
% qcond=0.58;
% hcond=0.6;
% regime='trans';



% Function handles

% B=@(x) exp(-abs((10.*(x-0.5))).^2)+1.5;   %Batimetria
% B=@(x) (sin(pi.*x./10)).^2;
% B=@(x) 1.2+0.*x;
 % B=@(x) 0.1+0.*x;
% B1=@(x) 0.05.*sin(x-12.5).*exp(1-(x-12.5).^2);
% B=@(x) B1(x)-min(B1(x));
% B=@(x) -0.06.*x.^2+1;
% B=@(x) 2*(sin(pi.*x)).^2;
% B=@(x) 0.05.*sin(x-12.5).*exp(1-(x-12.5).^2)+0.3;
% B=@(x) 0.2+0.5.*exp(1-1./(1-((abs(x-10)./5).^2))).*(x<15).*(x>5)+0*x;
% B=@(x) 0.*x+(0.2-0.05.*(x-10).^2).*(x<=12).*(x>=8);


% h0=@(x,xmax,xmin) 1+0.1.*sin(2.*pi.*x./(xmax-xmin));
% h0=@(x) w(x)-B(x);
% h0=@(x) 2.8-B(x);
% h0=@(x) 5+exp(cos(2.*pi.*x));
% h0=@(x) (2.1.*(x>0)+1.1.*(x<=0))-B(x);
% h0=@(x) (3+0.*x+2.5.*(x<0.6).*(x>0.4));
% h0=@(x) 2.*(exp(-abs(80.*(x-0.5)).^9)+2)-B(x);
% h0=@(x) (3-B(x)).*(x<10)+(1.5-B(x)).*(x>=10);
% q0=@(x) -1.*(x<=0.5)+1.*(x>0.5);
% q0=@(x) sin(cos(2.*pi.*x));
q0=@(x) 0.*x;


%Per le condizioni al bordo
% regime='lar'; 
% regime='super';
% regime='sub';
% regime='trans';
% qcond=1;
% hcond=1;

w=@(x) h0(x)+B(x);                          % w=h+B

% Discretizzazione spaziale
Nx=204;                         % Numero di celle
xmin=0; xmax=25;                % Bordi del dominio
x=linspace(xmin,xmax,Nx);             % Posizioni dei centri delle celle
dx=diff(x); dx=dx(1);           % delta x

% Discretizzazione temporale
Nt=5000;                         % Numero di step in tempo


% Parametri fisici
g=9.80665;                      % Costante di accelerazione gravitazionale
n=0.0;                         % Coefficiente di attrito di manning

% Inizializzazione variabili
h=zeros(Nt,Nx);                 % Altezza della colonna d'acqua
q=zeros(Nt,Nx);                 % Portata = h*u, con u velocità media

% Assegno a tutta la cella i valori sul baricentro
h(1,:)=h0(x);
q(1,:)=q0(x);


%CORRISPETTIVI PER IL METODO LAX-FRIEDRICHS
d=h;
v=q;




switch regime
    case 'lar' 
    case 'super'
        q(1,1)=qcond;
        h(1,1)=hcond;

        v(1,1)=qcond;
        d(1,1)=hcond;
    case 'sub'
        q(1,1)=qcond;
        h(1,end)=hcond;

        v(1,1)=qcond;
        d(1,end)=hcond;
    case 'trans'
        q(1,1)=qcond;
        h(1,end)=min(h(1,end),hcond);

        v(1,1)=qcond;
        d(1,end)=min(h(1,end),hcond);
end


% Indicizzazione con condizioni di Neumann
I=1:Nx;
R=I+1; R(end)=I(end);
% R2=R+1; R2(end-1)=1;
L=I-1; L(1)=I(1);
% L2=L-1; L(2)=Nx;
% dt=dx;
% dt=dx*0.1;



% GRAFICA
y_acqua=h(1,:)+B(x);
y_base = zeros(size(x)); % La base dell'area è a y=0
y_patch = [y_base, fliplr(y_acqua)]; % Combina base e f2 per creare il patch; IN REALTA' PARTE DA 0 IL BLOCCO IN ALTEZZA

% Normalizza i valori di y per la scala di colori
colors = (y_patch - min(y_patch)) / (max(y_patch) - min(y_patch)); % Normalizzati tra 0 e 1


% % Crea un oggetto VideoWriter per salvare l'animazione in MP4
% v = VideoWriter('animazione.mp4', 'MPEG-4');  % Specifica il nome del file e il formato
% v.FrameRate = 60;  % Imposta il frame rate (esempio: 10 FPS)
% open(v);  % Apre il file video per scrittura


%FIGURA PER IL METODO CON GLOBAL FLUX
figure('MenuBar', 'none', 'ToolBar', 'none');


% c sere per non vedere parte degli assi
y_max_vect=(max(w(x))+4)+0.*x;
c = patch('XData', [x, fliplr(x)], ... % Vertici x (andata e ritorno)
          'YData', [y_base, y_max_vect] , ... % Vertici y
          'FaceVertexCData', colors(:), ... % Scala di colori verticali
          'FaceColor', [0.91, 0.96, 1], ... % Interpolazione dei colori
          'EdgeColor', 'none', ...
          'linewidth', eps); % Senza bord

hold on;

% Crea il patch per l'area con scala di colori verticale
p = patch('XData', [x, fliplr(x)], ... % Vertici x (andata e ritorno)
          'YData', y_patch, ... % Vertici y
          'FaceVertexCData', colors(:), ... % Scala di colori verticali
          'FaceColor', 'interp', ... % Interpolazione dei colori
          'EdgeColor', 'none', ...
          'linewidth', eps); % Senza bordi
% hold on;
% Area sotto f1(x) (marrone scuro)
l = area(x, B(x), 'FaceColor', [0.3, 0.2, 0], 'LineWidth', 2.5, 'EdgeColor', [0.28, 0.18, 0]);

% Personalizzazione del grafico
% axis([x(1) x(end) 0 max(max(f1(x)),max(f2(x))+2)]);
axis([x(1) x(end) 0 max(c.YData)]);
grid off;
set(gca, 'Color', [0.91, 0.96, 1]); % Sfondo azzurro chiaro
% Aggiungi colormap blu
colors=150;
map = [(linspace(0,0.76,colors))',(linspace(0,0.9,colors))',(linspace(0.8,1,colors))']; % I COLORI PIU' SCURI SONO COPERTI DAL FONDALE
colormap(map); % Scala verticale di blu

set(gcf, 'InvertHardcopy', 'on'); % Forza lo sfondo bianco durante l'esportazione
set(gcf, 'Color', 'w'); % Sfondo bianco della finestra figura

% Cattura il frame corrente
ax=gca;
    ti = ax.TightInset;  % Calcola i margini stretti degli assi
    set(ax, 'Position', [ti(1), ti(2), 1 - ti(1) - ti(3), 1 - ti(2) - ti(4)]);

    % Imposta le dimensioni della figura
    set(gcf, 'PaperPositionMode', 'auto'); % Usa la dimensione corrente
    set(gcf, 'InvertHardcopy', 'off');     % Mantieni i colori della figura
    % frame = getframe(gcf);  % Cattura la figura corrente
    % % Estrai l'immagine dal frame (è un campo chiamato 'cdata')
    % img = frame.cdata;
    % 
    % % Crea il nome del file per il frame corrente
    % filename = sprintf('frameLAR.png');
    % % Salva l'immagine come file PNG
    % imwrite(img, filename);



% FIGURA PER IL CONFRONTO TRA I METODI
figure('Name','2','MenuBar', 'none', 'ToolBar', 'none')
d_patch = [y_base, fliplr(y_acqua)]; % Combina base e f2 per creare il patch; IN REALTA' PARTE DA 0 IL BLOCCO IN ALTEZZA

% Normalizza i valori di y per la scala di colori
dcolors = (d_patch - min(d_patch)) / (max(d_patch) - min(d_patch)); % Normalizzati tra 0 e 1


% % Crea un oggetto VideoWriter per salvare l'animazione in MP4
% v = VideoWriter('animazione.mp4', 'MPEG-4');  % Specifica il nome del file e il formato
% v.FrameRate = 60;  % Imposta il frame rate (esempio: 10 FPS)
% open(v);  % Apre il file video per scrittura


% c sere per non vedere parte degli assi
y_max_vect=(max(w(x))+4)+0.*x;
cielo = patch('XData', [x, fliplr(x)], ... % Vertici x (andata e ritorno)
          'YData', [y_base, y_max_vect] , ... % Vertici y
          'FaceVertexCData', dcolors(:), ... % Scala di colori verticali
          'FaceColor', [0.91, 0.96, 1], ... % Interpolazione dei colori
          'EdgeColor', 'none', ...
          'linewidth', eps); % Senza bord

hold on;

% Crea il patch per l'area con scala di colori verticale
dplot = patch('XData', [x, fliplr(x)], ... % Vertici x (andata e ritorno)
          'YData', d_patch, ... % Vertici y
          'FaceVertexCData', dcolors(:), ... % Scala di colori verticali
          'FaceColor', 'interp', ... % Interpolazione dei colori
          'EdgeColor', 'none', ...
          'linewidth', eps); % Senza bordi
% hold on;
% Area sotto f1(x) (marrone scuro)
bati = area(x, B(x), 'FaceColor', [0.3, 0.2, 0], 'LineWidth', 2.5, 'EdgeColor', [0.28, 0.18, 0]);

vol=plot(x,h(1,:)+B(x),'r');

% Personalizzazione del grafico
% axis([x(1) x(end) 0 max(max(f1(x)),max(f2(x))+2)]);
axis([x(1) x(end) 0 max(c.YData)]);
grid off;
ax=gca;
    ti = ax.TightInset;  % Calcola i margini stretti degli assi
    set(ax, 'Position', [ti(1), ti(2), 1 - ti(1) - ti(3), 1 - ti(2) - ti(4)]);

    % Imposta le dimensioni della figura
    set(gcf, 'PaperPositionMode', 'auto'); % Usa la dimensione corrente
    set(gcf, 'InvertHardcopy', 'off');     % Mantieni i colori della figura
set(gca, 'Color', [0.91, 0.96, 1]); % Sfondo azzurro chiaro
% Aggiungi colormap blu
colors=150;
map = [(linspace(0,0.76,colors))',(linspace(0,0.9,colors))',(linspace(0.8,1,colors))']; % I COLORI PIU' SCURI SONO COPERTI DAL FONDALE
colormap(map); % Scala verticale di blu

set(gcf, 'InvertHardcopy', 'on'); % Forza lo sfondo bianco durante l'esportazione
set(gcf, 'Color', 'w'); % Sfondo bianco della finestra figura



% FIGURA PER IL CONFRONTO DELLE q
figure('Name','3','MenuBar', 'none', 'ToolBar', 'none')
q_fv=plot(x,q(1,:),'r',LineWidth=1.5);
hold on;
q_lf=plot(x,v(1,:),'b',LineWidth=1.5);
serve_per_tempo=plot(nan,nan,'w');
axis([x(1) x(end) -0.2 0.2]);
grid off;
ax=gca;
    ti = ax.TightInset;  % Calcola i margini stretti degli assi
    set(ax, 'Position', [ti(1), ti(2), 1 - ti(1) - ti(3), 1 - ti(2) - ti(4)]);

    % Imposta le dimensioni della figura
    set(gcf, 'PaperPositionMode', 'auto'); % Usa la dimensione corrente
    set(gcf, 'InvertHardcopy', 'off');     % Mantieni i colori della figura
set(gcf, 'InvertHardcopy', 'on'); % Forza lo sfondo bianco durante l'esportazione
set(gcf, 'Color', 'w'); % Sfondo bianco della finestra figura


t=0;

tempo=sprintf('t= %.2f',t);
leg=legend({'q con volumi finiti','q con Lax-Friedrichs',tempo});

      
%FIGURA PER IL CONFRONTO DELLE h
figure('Name','4','MenuBar', 'none', 'ToolBar', 'none')
h_fv=plot(x,h(1,:)+B(x),'r',LineWidth=1.5);
hold on;
h_lf=plot(x,d(1,:)+B(x),'b',LineWidth=1.5);
serve_per_tempo=plot(nan,nan,'w');
axis([x(1) x(end) 2.14 2.26]);
grid off;
ax=gca;
    ti = ax.TightInset;  % Calcola i margini stretti degli assi
    set(ax, 'Position', [ti(1), ti(2), 1 - ti(1) - ti(3), 1 - ti(2) - ti(4)]);

    % Imposta le dimensioni della figura
    set(gcf, 'PaperPositionMode', 'auto'); % Usa la dimensione corrente
    set(gcf, 'InvertHardcopy', 'off');     % Mantieni i colori della figura
set(gcf, 'InvertHardcopy', 'on'); % Forza lo sfondo bianco durante l'esportazione
set(gcf, 'Color', 'w'); % Sfondo bianco della finestra figura


t=0;

tempo=sprintf('t= %.2f',t);
leg1=legend({'h con volumi finiti','h con Lax-Friedrichs',tempo});

                
waitforbuttonpress;
eta=h(1,:)+B(x);

tol=10^(-6);
k=1;
norma=tol+1; %per entrare nel primo while
perturbato=0; %non è stata introdotta alcuna perturbazione

%PER PERTURBARE UN LAR COMMENTARE LA SECONDA CONDIZIONE DEL WHILE
while (k<Nt) & (norma>tol)
     % velocità media
    u=q(k,:)./h(k,:);
    
    % CFL
    c=max(abs(u)+(g.*h(k,:)).^(1/2));
    dt=dx/c*0.5;
    t=t+dt;
        tempo=sprintf('t= %.2f',t);
        leg.String = {'q con volumi finiti','q con Lax-Friedrichs',tempo};
        
        leg1.String = {'h con volumi finiti','h con Lax-Friedrichs',tempo};

    %trovo K1

    h_old=h(k,:);
    q_old=q(k,:);
    % B_app=(B(x)+B(x(L)))./2;
    % B_app2=(B(x)+B(x(R)))./2;
    [Kh1,Kq1]=rk_step(h_old,q_old,Nx,dx,L,R,g,B,x,n);


    % Trovo K2

    h_old=h(k,:)+dt.*0.5.*Kh1;
    q_old=q(k,:)+dt.*0.5.*Kq1;
    [Kh2,Kq2]=rk_step(h_old,q_old,Nx,dx,L,R,g,B,x,n);
    


     % Trovo K3

    h_old=h(k,:)+dt.*0.5.*Kh2;
    q_old=q(k,:)+dt.*0.5.*Kq2;
    [Kh3,Kq3]=rk_step(h_old,q_old,Nx,dx,L,R,g,B,x,n);
    


     % Trovo K4

    h_old=h(k,:)+dt.*Kh3;
    q_old=q(k,:)+dt.*Kq3;
    [Kh4,Kq4]=rk_step(h_old,q_old,Nx,dx,L,R,g,B,x,n);
    


    h(k+1,:)=h(k,:)+dt.*(Kh1+Kh2+Kh3+Kh4)./6;
    q(k+1,:)=q(k,:)+dt.*(Kq1+Kq2+Kq3+Kq4)./6;

    
    %discretizzazione della sorgente per L-F
    S=-0.5.*g.*d(k,:).*(B(x(R))-B(x(L)))./dx;
    for i=1:Nx
        % vel=v(k,i)./d(k,i);
        % if vel>=0
        %     d(k+1,i)=d(k,i)-dt./dx.*(v(k,i)-v(k,L(i)));
        %     v(k+1,i)=v(k,i)-dt.*(v(k,i).^2./d(k,i)+g.*0.5.*d(k,i).^2-v(k,L(i)).^2./d(k,L(i))-g.*0.5.*d(k,L(i)).^2)+dt.*S(i);
        % elseif vel<0
        %     d(k+1,i)=d(k,i)-dt./dx.*(-v(k,i)+v(k,R(i)));
        %     v(k+1,i)=v(k,i)-dt.*(-v(k,i).^2./d(k,i)-g.*0.5.*d(k,i).^2+v(k,R(i)).^2./d(k,R(i))+g.*0.5.*d(k,R(i)).^2)+dt.*S(i);
        % end
        d(k+1,i)=0.5.*(d(k,R(i))+d(k,L(i)))-0.5.*dt./dx.*(v(k,R(i))-v(k,L(i)));
        v(k+1,i)=0.5.*(v(k,R(i))+v(k,L(i)))-0.5.*dt./dx.*(v(k,R(i)).^2./d(k,R(i))+0.5.*g.*d(k,R(i)).^2-v(k,L(i)).^2./d(k,L(i))-0.5.*g.*d(k,L(i)).^2)+dt.*S(i);
    end


    switch regime
        case 'lar' 
        case 'super'
            q(k+1,1)=qcond;
            h(k+1,1)=hcond;

            v(k+1,1)=qcond;
            d(k+1,1)=hcond;
        case 'sub'
            q(k+1,1)=qcond;
            h(k+1,end)=hcond;


            v(k+1,1)=qcond;
            d(k+1,end)=hcond;
        case 'trans'
            q(k+1,1)=qcond;
            h(k+1,end)=min(h(k+1,end),hcond);


            v(k+1,1)=qcond;
            d(k+1,end)=min(d(k+1,end),hcond);
    end
    
      
    %per la grafica
    y_patch = [0.*x, fliplr(h(k+1,:)+B(x))]; % PER FARLO PARTIRE DAL FONDALE E NON DA 0 CAMBIARE LA PRIMA COMPONENTE IN q.YData
    % Normalizza i valori di y per la scala di colori
    colors = (y_patch - min(y_patch)) / (max(y_patch) - min(y_patch)); % Normalizzati tra 0 e 1

    p.YData=y_patch;
    p.FaceVertexCData=colors(:);
    % c.YData= [max(y_patch(end/2+1:end),l.YData), (max(w(x))+2)+0.*x];
    % pausa=0; %serve a personalizzare la pausa tra i frame prima e dopo le perturbazioni
    if perturbato==1
        oscill.YData=h(k+1,:)+B(x);
        oscill2.YData=d(k+1,:)+B(x);
        % pausa=0;
    end
    
    

    %per la grafica
    d_patch = [0.*x, fliplr(d(k+1,:)+B(x))]; 
    % Normalizza i valori di y per la scala di colori
    dcolors = (d_patch - min(d_patch)) / (max(d_patch) - min(d_patch)); % Normalizzati tra 0 e 1
    dplot.YData=d_patch;
    dplot.FaceVertexCData=dcolors(:);
    vol.YData=h(k+1,:)+B(x);

    q_fv.YData=q(k+1,:);
    q_lf.YData=v(k+1,:);

    h_fv.YData=h(k+1,:)+B(x);
    h_lf.YData=d(k+1,:)+B(x);

    drawnow;
    % pause(pausa);

    % % Cattura il frame corrente
    % frame = getframe(gcf);  % Cattura la figura corrente
    % 
    % writeVideo(v, frame);  % Scrivi il frame nel video
    
    %per il criterio d'arresto
    norma=max(abs(h(k+1,:)-h(k,:)));

    %PER LE PERTURBAZIONI; SCOMMENTARE QUESTA SEZIONE: PARTE UNA VOLTA RAGGIUNTO UNO STATO DI QUASI-EQUILIBRIO DATO DALLA CONDIZIONE DEL WHILE
    % if (perturbato==0) & (norma<10.*tol)
    %     figure('Name', '3','MenuBar', 'none', 'ToolBar', 'none');
    %     equilibrio=plot(x,h(k+1,:)+B(x),'g', LineWidth=1);
    %     hold on;
    %     pert=0.05.*(x>3.8).*(x<4.5);
    %     h(k+1,:)=h(k+1,:)+pert;
    %     d(k+1,:)=d(k+1,:)+pert;
    %     oscill=plot(x,h(k+1,:)+B(x),'r',LineWidth=1);
    %     oscill2=plot(x,d(k+1,:)+B(x),'b',LineWidth=1);
    %     leg=legend('Equilibrio','Dinamica perturbata','Dinamica perturbataLax-Friedrichs');
    %     axis([x(1) x(end) 2.18 2.3]);
    %     set(gcf, 'InvertHardcopy', 'on'); % Forza lo sfondo bianco durante l'esportazione
    %     set(gcf, 'Color', 'w'); % Sfondo bianco della finestra figura
    %     leg.Position=[0.62,0.4, 0.25 ,0.05];
    %     y_patch = [0.*x, fliplr(h(k+1,:)+B(x))]; % PER FARLO PARTIRE DAL FONDALE E NON DA 0 CAMBIARE LA PRIMA COMPONENTE IN q.YData
    %     % Normalizza i valori di y per la scala di colori
    %     colors = (y_patch - min(y_patch)) / (max(y_patch) - min(y_patch)); % Normalizzati tra 0 e 1
    %     d_patch = [0.*x, fliplr(d(k+1,:)+B(x))]; % PER FARLO PARTIRE DAL FONDALE E NON DA 0 CAMBIARE LA PRIMA COMPONENTE IN q.YData
    %     % Normalizza i valori di y per la scala di colori
    %     dcolors = (d_patch - min(d_patch)) / (max(d_patch) - min(d_patch)); % Normalizzati tra 0 e 1
    % 
    %     p.YData=y_patch;
    %     p.FaceVertexCData=colors(:);
    %     dplot.YData=d_patch;
    %     dplot.FaceVertexCData=dcolors(:);
    %     drawnow;
    %     pause(4);
    %     perturbato=1;
    % end
    k=k+1;
end

% Chiudi il file video
% close(v);



function [Kh,Kq]=rk_step(h,q,Nx,dx,L,R,g,B,x,n)
    JF=zeros(2,2,Nx);
    Lambda=zeros(2,2,Nx);
    Ri=zeros(2,2,Nx);
    Ri_inv=zeros(2,2,Nx);
    r=zeros(1,Nx);
    eta=h+B(x);
    % eta_2=B_app2+h_2;
    % eta_=B_app+h_;
    % eta_=(eta+eta(L))./2;
    for i=1:Nx
        if i>1
            r(i)=r(i-1)+g.*0.5.*(B(x(R(i-1)))-B(x(i-1))).*(eta(R(i-1))+eta(i-1))-0.5.*g.*(B(x(R(i-1))).^2-B(x(i-1)).^2);
            r(i)=r(i)+0.5.*dx.*n.*n.*(abs(q(i-1)).*q(i-1)./(h(i-1).^(7/3))+abs(q(i)).*q(i)./(h(i).^(7/3))); %aggiunta parametro con la frizione
            % r(i)=r(i-1)+g.*0.5.*(h(i)+h(i-1)).*(B(x(i))-B(x(i-1))); %altra formula per il termine senza frizione
        end
        h_star=(h(i)+h(R(i)))/2;
        % q_star=(q(i)+q(R(i)))/2;    % questa q_star aumenta gli errori
        q_star=((h(i).^0.5.*q(i)./h(i)+h(R(i)).^0.5.*q(R(i))./h(R(i)))./(h(i).^0.5+h(R(i)).^0.5)).*h_star;
        ul=q(i)./h(i);
        ur=q(R(i))./h(R(i));
        % u_star=(sqrt(h(i)).*ul+sqrt(h(R(i))).*ur)./(sqrt(h(i))+sqrt(h(R(i))));
        
        %conto fatto con le funzioni matlab
        % JF(:,:,i)=[0,1;-q_star.^2./h_star.^2+g.*h_star, 2.*q_star./h_star];             %QUESTO CONTO VA FATTO AD OGNI STEP DI RK4?
        % [Ri(:,:,i),Lambda(:,:,i)]=eig(JF(:,:,i));
        % Ri_inv(:,:,i)=inv(Ri(:,:,i));

        %conto fatto a mano
        % Questa scelta è più rapida computazionalmente
        a=0;
        b=1;
        c=-q_star.^2./h_star.^2+g.*h_star;
        d=2.*q_star./h_star;
        JF(:,:,i)=[a,b;c,d];
        l1=(d-(d.^2+4.*c).^0.5)./2;
        l2=(d+(d.^2+4.*c).^0.5)./2;
        Lambda(:,:,i)=[l1,0;0,l2];
        Ri(:,:,i)=[1,1;l1,l2];
        Ri_inv(:,:,i)=1./(l2-l1).*[l2,-1;-l1,1];
        
    end
    
    % Calcolo della componente del campo relativa a q
    G=q.^2./(h)+0.5.*g.*h.^2+r;
    
    
    app=[q;G];
    appdx=[q(R);G(R)];

    for i=1:Nx
        v(:,i)=Ri(:,:,i)*(Lambda(:,:,i)>0)*Ri_inv(:,:,i)*app(:,i)+Ri(:,:,i)*(Lambda(:,:,i)<0)*Ri_inv(:,:,i)*appdx(:,i);
    end


    q_cap=v(1,:);
    G_cap=v(2,:);


    Kh=-(q_cap-q_cap(L))./dx;
    Kq=-(G_cap-G_cap(L))./dx;
end




