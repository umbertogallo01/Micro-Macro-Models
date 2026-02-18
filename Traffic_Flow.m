%% ES 1 - TRAFFIC FLOW MODELS (Follow-the-leader + Optimal velocity)
% A. Write a function for the numerical solution of the traffic flow model FTL + FOV 
% with both first– and second–order explicit schemes (forward Euler and, e.g., Heun). 
% The input arguments are:
% - the number of vehicles N
% - the initial conditions (x0_i, v0_i)
% - the length of the road L
% - the parameter α 
% - the parameter β 
% - the final observation time T
% - the time step ∆t
% - the order of the method
% 
% The output arguments are:
% - the position in time of each vehicle
% - the velocity in time of each vehicle

clc
clear
close all


% Parametri fissati per la simulazione numerica
N = 90;
L = 1500;
dt = 0.1;
T = 1000;


% Ordine dello schema numerico usato per il sistema di ODE
% ord = 1;   % Eulero Esplicito
ord = 2;    % Heun


% Parametri fissati del modello del traffico
alfa = 1;       
beta = 100;     
V1 = 6.75; 
V2 = 7.91;
C1 = 0.13;      % parametro di calibrazione 1
C2 = 1.57;      % parametro di calibrazione 2
l = 5;          % lunghezza dei veicoli
dx=L/N;         % distanza tra i veicoli al tempo iniziale (equispaziati)


% Condizioni iniziali
X0 = (dx:dx:L)';
V0 = Velocity(dx.*ones(N,1),V1,V2,C1,C2,l);


% Modello Follow-the-leader + Bando
[x,v,pho,t] = traffic_flow_model(N,X0,V0,L,alfa,beta,T,dt,ord,V1,V2,C1,C2,l);


% Numero di passi della discretizzazione
Nt=length(t);


% Visualizzazione grafica dei risultati: supponiamo di rappresentare tutto
% o su una strada circolare lunga L, quindi trasformiamo le posizioni dei 
% veicoli in coordinate polari o su una strada rettilinea

strada=menu('Scegliere il tipo di strada: ', 'Circolare', 'Rettilinea');

switch strada
    case 1
        theta = (2*pi/L).*x;
        r = L/(2*pi);
        X1 = r.*cos(theta);
        X2 = r.*sin(theta);

    case 2
        x=mod(x,L); % Per avere condizioni al bordo periodiche
end


% 1) Plot della velocità del veicolo 1 e delle traiettorie (t,x) di tutti i
% veicoli

figure(1)
set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero

subplot(1,2,1)
p1 = plot(t,v(1,:)); 
xlabel 'Tempo'
ylabel 'Velocità'
grid on
%ylim([0 14.66]); % V=14.66 m/s è la massima velocità
ylim([min(v,[],"all")-1 max(v,[],"all")+1])
title('Evoluzione della velocità del veicolo 1','FontSize',14);
axis square

subplot(1,2,2)
f1=plot(t,x);
xlabel 'Tempo';
ylabel 'Posizioni';
xlim([0 T])
title('Traiettorie di tutti i veicoli','FontSize',14);
grid on
axis square;


% 2) I veicoli sono inizialmente equidistribuiti e supponiamo che 
% v0_i = V (L/N), N=120, L=1500 e T=1000

% Test 2a: siano α = 1, β = 0 (Bando Model). Verificare che le condizioni 
% iniziali definiscono uno stato di equilibrio del modello. 
% Successivamente supporre di introdurre una perturbazione al tempo t=100 
% del veicolo i=N/2 che riduce la sua velocità a v_N/2(100)=ν*V(L/N), ν∈(0,1];

% Test 2b: rifare la stessa cosa precedente con β = 100

waitforbuttonpress;
figure(2)
set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero

switch strada

    case 1
        ax_anim = gca; % ottieni il riferimento all'asse corrente
        f2 = plot(X1(:,1), X2(:,1),'o','MarkerSize',5,'MarkerFaceColor','#D9FFFF');
        
        hold on
            f_veicolo_rall = plot(X1(N/2,1),X2(N/2,1),'o','MarkerSize',5,'MarkerFaceColor','g');
            f_veicolo_autonomo = plot(X1(1,1),X2(1,1),'o','MarkerSize',5,'MarkerFaceColor','r');
        hold off

        xlabel('X1')
        ylabel('X2')
        %legend('Veicolo rosso = veicolo autonomo', 'Veicolo verde = veicolo che frena')
        title(ax_anim, ['Tempo = ', num2str(0)]);

        grid on
        axis square; 

        % Ciclo per l'animazione
        for n = 1:3:Nt-1
    
            f2.XData = X1(:,n+1);
            f2.YData = X2(:,n+1);   
        
            f_veicolo_rall.XData = X1(N/2,n+1);
            f_veicolo_rall.YData = X2(N/2,n+1);
        
            f_veicolo_autonomo.XData = X1(1,n+1);
            f_veicolo_autonomo.YData = X2(1,n+1);
        
        
            % Rende visibile l'aggiornamento
            drawnow;
        
            % Aggiorna solo il titolo dell’asse della figura animata
            title(ax_anim, ['Tempo = ', num2str(t(n))]);

            pause(0.01);

        end

    case 2
        ax_anim = gca; % ottieni il riferimento all'asse corrente
        y= 500*ones(N,1);
        f2 = plot(x(:,1), y  ,'o','MarkerSize',5,'MarkerFaceColor','#D9FFFF');

        yline(500.05, 'k', 'LineWidth', 2);  % linea orizzontale a y = 3
        yline(499.95, 'k', 'LineWidth', 2);  % linea orizzontale a y = 7
        
        hold on
            f_veicolo_rall = plot(x(N/2,1),500,'o','MarkerSize',5,'MarkerFaceColor','g');
            f_veicolo_autonomo = plot(x(1,1),500,'o','MarkerSize',5,'MarkerFaceColor','r');
        hold off

        xlabel('Strada')
        ylabel('Veicoli')
       % legend('Veicolo rosso = veicolo autonomo', 'Veicolo verde = veicolo che frena')
        ylim([499.95-1 500.05+1])
        title(ax_anim, ['Tempo = ', num2str(0)]);

        % Ciclo per l'animazione
        for n = 1:5:Nt-1

            f2.XData = x(:,n+1);   

            f_veicolo_rall.XData = x(N/2,n+1);
        
            f_veicolo_autonomo.XData = x(1,n+1);
        
        
            % Rende visibile l'aggiornamento
            drawnow;
        
            % Aggiorna solo il titolo dell’asse della figura animata
            title(ax_anim, ['Tempo = ', num2str(t(n))]);

            pause(0.01);

        end

end


% 3) I veicoli sono inizialmente equidistribuiti e supponiamo che 
% v0_i = V (L/N), N=90, L=1500 e T=1000. Ripetere i test del punto 2)


% 4) Mostrare l'evoluzione nel tempo della densità locale e della velocità
% di ogni veicolo

figure(3)
subplot(1,2,1)
ax_anim1 = gca; % ottieni il riferimento all'asse corrente

density=plot(X0,pho(:,1));

set(gcf, 'WindowState', 'maximized');   % Imposta la figura a schermo intero
xlabel('Strada')
ylabel('ρ(t)')
ylim([-0.2,1.2]);
title(ax_anim1, ['Evoluzione nel tempo della densità delle auto sulla strada (t = ' num2str(0) , ')']);
grid on
axis square;


subplot(1,2,2)
ax_anim2 = gca; % ottieni il riferimento all'asse corrente

velocity=plot(t(1),v(:,1));
set(gcf, 'WindowState', 'maximized');   % Imposta la figura a schermo intero
xlabel('Tempo');
ylabel('Velocità');
%xlim([0 L]);
ylim([V1-V2 V1+V2])
title(ax_anim2, ['Evoluzione nel tempo della velocità dei veicoli (t = ' num2str(0) , ')']);
grid on
axis square;


% Ciclo per l'animazione
for n = 1:10:Nt

    density.YData=pho(:,n+1);
    title(ax_anim1, ['Evoluzione nel tempo della densità delle auto sulla strada (t = ' num2str(t(n)) , ')']);

    for i=1:N
        velocity(i).XData=[velocity(i).XData t(n+1)];
        velocity(i).YData=[velocity(i).YData v(i,n+1)];
    end
   
    xlim(ax_anim2,[t(n)-20, t(n)])
    title(ax_anim2, ['Evoluzione nel tempo della velocità dei veicoli (t = ' num2str(t(n)) , ')']); %modifico gli assi solo del secondo subplot della figura 2

    % Rende visibile l'aggiornamento
    drawnow;

    pause(0.01);
end


%% Funzioni usate

function [x,v,pho,t] = traffic_flow_model(N,X0,V0,L,alfa,beta,T,dt,ord,V1,V2,C1,C2,l)
    
    % Controllo sul passo di discretizzazione in tempo
    t = 0:dt:T;

    if t(end) < T           % non ricopro [0,T] con il dt scelto
        Nt = length(t)+1;
    
    else                    % ricopro [0,T] con dt
        Nt = length(t);
    end

    dt = T/(Nt-1);
    t=0:dt:T;
  
    
    % Salviamo in un vettore y=(x,v) la soluzione per colonne
    y0 = [X0;V0];
    y = zeros(2*N,Nt);
    y(:,1) = y0;
    

    % Inizializziamo un vettore per le densità locali che avrà per ogni
    % colonna la densità rho_i al tempo t. Nella prima colonna, dato che
    % stiamo assumendo che tutti i veicoli siano equispaziati, avremmo lo
    % stesso headway

    c = [y(2:N,1); L + y(1,1)]; % periodic boundary conditions
    d = y(1:N,1); 
    dx = c-d;
    pho=1./dx;
    

    n=1;
    while t(n)<T

        % PERTURBAZIONE: al tempo t=100, supponiamo che il veicolo i=N/2
        % freni riducendo la sua velocità a v_N/2(100)=ν * V(L/N), ν∈(0,1];

        if t(n)==100

            nu=0.001;
            y(N+(N/2),n) = nu * Velocity(L/N,V1,V2,C1,C2,l);

        end


        % Risolvo il sistema di ODE : y'(t)=F(t,y(t))
        switch ord 
            
            case 1 % Eulero Esplicito

                [Fn,dx] = F(y(:,n), N, alfa, beta, t(n), L, V1, V2, C1, C2, l);
                y(:,n+1) = y(:,n) + dt*Fn;  % Eulero esplicito
                pho(:,n+1)=l./dx;
                n=n+1;

            case 2 % Heun

               [Fn1,dx] = F(y(:,n), N, alfa, beta, t(n), L, V1, V2, C1, C2, l);
               Fn2 = F(y(:,n)+dt*Fn1, N, alfa, beta, t(n), L, V1, V2, C1, C2, l);
               y(:,n+1) = y(:,n) + (dt/2)*(Fn1 + Fn2);  % Heun
               pho(:,n+1)=l./dx;
               n=n+1;

        end

    end

    % Estraggo le posizioni: y(:,t)=[x1(t),...,xN(t),v1(t),...,vN(t)]
    x = y(1:N,:);

    % Estraggo le velocità: y(:,t)=[x1(t),...,xN(t),v1(t),...,vN(t)] 
    v = y(N+1:2*N,:); 
end



% Funzione di velocità: la velocità minima attesa è data da V1-V2, mentre
% quella massima V1+V2. Dunque nel nostro caso la velocità massima 
% è Vmax=14.66 m/s
function V = Velocity(dx,V1,V2,C1,C2,l)

    r = length(dx);                         % numero di righe di dx
    V_tilde = V1 + V2*tanh(C1*(dx-l)-C2);
    V = max(V_tilde, zeros(r,1));

end



% Dinamica del modello Follow-The-Leader + Bando: stiamo pensando al
% sistema in forma vettorizzata y'(t)=F(t,y(t)), dunque 
% F=[v ; beta*(Dv/Dx^2) + alpha * (V(Dx)-v)]

function [fx,dx] = F(y,N,alfa,beta,t,L,V1,V2,C1,C2,l)
    
    % Usiamo condizioni al bordo periodiche. Supponiamo quindi  di
    % considerare un circuito chiuso di lunghezza L e che la
    % la velocità del veicolo davanti a N sia uguale a quella del 
    % veicolo 1, mentre la sua posizione è data da L+x1(t)
    

    % VEICOLO AUTONOMO: supponiamo che al tempo t>=300, il veicolo 1 inizi 
    % a muoversi in modo controllato con una dinamica diversa. (Nota: per
    % studiare cosa succede senza veicolo autonomo, basta mettere M=1)

    if t>=300
        
        % Dinamica veicolo autonomo

        M=60; % Numero di veicolo da cui dipende il veicolo autonomo

        % calcoliamo dv, dx e v e V(dx) per scrivere la dinamica

        a = [y(N+2:end) ; y(N+1)];    % a = [v2(tn),...,vN(tn),v1(tn)]
        b = y(N+1:end);               % b = [v1(tn),...,vN-1(tn),vN(tn)]
        dv = a-b;
        

        c = [y(2:N) ; L + y(1)];      % Condizioni al bordo periodiche
        d = y(1:N); 
        dx = c-d;
       
        v = y(N+1:end);
        V = Velocity(dx(2:N),V1,V2,C1,C2,l);

        % Campo vettoriale del modello del traffico con veicolo autonomo
        fx = [v; beta*(dv(1)/(dx(1)^2)) + alfa*( Velocity((1/M)*sum(dx(1:M)),V1,V2,C1,C2,l) - v(1)) ;  beta*( dv(2:N)./((dx(2:N)).^2) ) + alfa*(V - v(2:N)) ];

    else
        
        % calcoliamo dv, dx e v e V(dx) per scrivere la dinamica

        a = [y(N+2:end) ; y(N+1)];    % a = [v2(tn),...,vN(tn),v1(tn)]
        b = y(N+1:end);               % b = [v1(tn),...,vN-1(tn),vN(tn)]
        dv = a-b;
        

        c = [y(2:N) ; L + y(1)];      % Condizioni al bordo periodiche
        d = y(1:N); 
        dx = c-d;
    
        v = y(N+1:end);
        V = Velocity(dx,V1,V2,C1,C2,l);
        

        % Campo vettoriale del modello del traffico senza veicolo autonomo
        fx = [v ; beta*(dv./((dx).^2)) + alfa*(V-v)];

    end

end
