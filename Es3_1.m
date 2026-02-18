%% ES 1 Foglio 3- WAVE EQUATION
% Write a code which numerically solves (1), endowed with an initial 
% condition U(x, t = 0) = U0(x), with an explicit first order finite volume
% scheme. In particular, consider the explicit Euler method for the time 
% discretization combined with the Rusanov numerical flux computed on
% piecewise constant reconstruction of the solution.

clc
clear
close all

disp("*** WAVE EQUATION ***")


% Parametri fissati per lo schema numerico
c=1;
A=[0 c; c 0];
xm=-1;
xM=1;
T=2;
N=200;                  % numero di celle
dx=(xM-xm)/N;
x=xm+((1:N)-0.5).*dx;   % centri delle celle


% Calcolo del time-step con le condizioni di stabilità CFL
lambda=max(abs(eig(A)));
C=0.5;  % costante di sicurezza
dt=(dx/lambda) * C;

t = 0:dt:T;
if t(end) < T           % non ricopro [0,Tfin] con dt della CFL 
    Nt = length(t)+1;

else                    % ricopro [0,Tfin] con dt della CFL
    Nt = length(t);
end
dt = T/(Nt-1);
t=0:dt:T;


% Inizializzazione dell'approssimazione numerica
scelta = menu('Scegli il dato iniziale', '1', '2', '3');
U=initial_conditions(x,N,Nt,scelta);


% Schema a volume finito
n=1;

while t(n)<T

    % Calcolo dei flussi di Rusanov
    Fd=Rusanov_flux(U(5:2*(N+2),n), U(3:2*(N+1),n), lambda, A, N);
    Fs=Rusanov_flux(U(3:2*(N+1),n), U(1:2*N,n), lambda, A, N);
    
    % Schema di Eulero esplicito
    U(3:2*(N+1),n+1) = U(3:2*(N+1),n) - (dt/dx)*(Fd-Fs);
    
    if scelta==1 || scelta==2      
        % Riapplichiamo le condizioni al bordo periodiche
        
        U(1:2,n+1)=U(end-3:end-2,n+1);   % Cella ghost di sinistra
        U(end-1:end,n+1)=U(3:4,n+1);     % Cella ghost di destra
        
    else
        % Riapplichiamo le condizioni al bordo di free-flow

        U(1:2,n+1)=U(3:4,n+1);                  % Cella ghost di sinistra
        U(end-1:end,n+1)=U(end-3:end-2,n+1);    % Cella ghost di destra
    end
      
    n=n+1;
end


% Visualizzazione grafica dei risultati

% Animazione
figure(1)
f1=plot(x, U(3:2:2*N+2,1), 'r', x, U(4:2:2*N+2,1), 'g');
set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero
grid on
axis square
xlabel('x')
ylabel('U(x,t)')
legend('v(x,t)','w(x,t)')
title(['Tempo = ', num2str(0)]);
ylim([min(U,[],"all")-1 max(U,[],"all")+1])

waitforbuttonpress
for i = 1:n-1

    f1(1).YData=U(3:2:2*N+2,i+1);
    f1(2).YData=U(4:2:2*N+2,i+1);

    % Rende visibile l'aggiornamento
    drawnow;

    % aggiorna solo il titolo dell’asse della figura animata
    title(['Tempo = ', num2str(t(i+1))]);
    pause(0.001);
end



figure(2)
if scelta==1 || scelta==2
    subplot(3,1,1);
    i=find(t==0.5);     % Individuo l'indice del vettore dei tempi corrispondente a t=0.5
    plot(x,U(3:2:2*N+2,i),'r',x,U(4:2:2*N+2,i),'g');
    set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero
    grid on
    xlabel('x')
    ylabel('U(x,t=0.5)')
    legend('v(x,t)','w(x,t)')
    title(['Tempo = ', num2str(0.5)]);

    subplot(3,1,2)
    i=find(t==1);       % Individuo l'indice del vettore dei tempi corrispondente a t=1
    plot(x,U(3:2:2*N+2,i),'r',x,U(4:2:2*N+2,i),'g');
    grid on
    xlabel('x')
    ylabel('U(x,t=1)')
    legend('v(x,t)','w(x,t)')
    title(['Tempo = ', num2str(1)]);

    subplot(3,1,3)
    i=find(t==2);       % Individuo l'indice del vettore dei tempi corrispondente a t=2
    plot(x,U(3:2:2*N+2,i),'r',x,U(4:2:2*N+2,i),'g');
    grid on
    xlabel('x')
    ylabel('U(x,t=2)')
    legend('v(x,t)','w(x,t)')
    title(['Tempo = ', num2str(2)]);
        
else
    subplot(2,1,1);
    i=find(t==0.5);     % Individuo l'indice del vettore dei tempi corrispondente a t=0.5
    plot(x,U(3:2:2*N+2,i),'r',x,U(4:2:2*N+2,i),'g');
    set(gcf, 'WindowState', 'maximized');
    grid on
    xlabel('x')
    ylabel('U(x,t=0.5)')
    legend('v(x,t)','w(x,t)')
    title(['Soluzione numerica per t = ', num2str(0.5)]);
    ylim([min(U,[],"all")-1 max(U,[],"all")+1])
    

    % Calcolo della soluzione esatta
    subplot(2,1,2)

    % Soluzione in variabili conservative tra x=-t e x=t
    R=ones(2);
    R(1,1)=-1;
    v=R*[(1-0)/2;(1+0.5)/2];
    
    % Definizione della soluzione esatta 
    V=@(x,t) 1.*(x<=-t) + v(1) .* (x>-t & x<t) + 0.*(x>=t);
    W=@(x,t) 0.5.*(x<=-t) + v(2) .* (x>-t & x<t) + 1.*(x>=t);

    % Grafico della soluzione esatta per t=0.5
    plot(x,V(x,0.5),'r',x,W(x,0.5),'g');
    
    grid on
    xlabel('x')
    ylabel('U esatta(x,t=0.5)')
    legend('v(x,t)','w(x,t)')
    title(['Soluzione esatta in variabili conservative per t = ', num2str(0.5)]);
    ylim([min(U,[],"all")-1 max(U,[],"all")+1])

     
    % % Soluzione esatta in variabili caratteristiche
    % figure(3)
    % P1=@(x,t) -0.5*V(x,t)+0.5*W(x,t);
    % P2=@(x,t) 0.5*V(x,t)+0.5*W(x,t);
    % plot(x,P1(x,0.5),'r',x,P2(x,0.5),'g');
    % set(gcf, 'WindowState', 'maximized');
    % xlabel('x')
    % ylabel('U_esatta(x,t=0.5)')
    % legend('v(x,t)','w(x,t)')
    % axis square
    % title(['Soluzione esatta in variabili caratteristiche per t = ', num2str(0.5)]);

end
        

%% Funzioni

% Condizioni iniziali: organizziamo l'approssimazione numerica come
% U=[Ghost_L , u0(x1),v0(x1), ... ,u0(xn),v0(xN) , Ghost_R]

function U = initial_conditions(x,N,Nt,scelta)
    
    % Inizializzazione dell'approssimazione numerica con le due celle ghost
    U=zeros(2*(N+2),Nt);
    
    switch scelta
        case 1
            u0=@(x) (1-(cos(2*pi*(x-0.25))).^2).*(x>=0.25 & x<=0.75); 

            % Definizione delle approssimazioni iniziali date dal profilo
            % iniziale valutato nei centri delle celle
            U(3:2:2*(N+1),1)=u0(x);
            
            % Creazione delle celle ghost con condizioni al bordo
            % periodiche
            U(1:2,1)=U(end-3:end-2,1);  % Cella ghost di sinistra
            U(end-1:end,1)=U(3:4,1);    % Cella ghost di destra


        case 2
            u0=@(x) 1.*(x>=0.25 & x<=0.75);
            
            % Definizione delle approssimazioni iniziali date dal profilo
            % iniziale valutato nei centri delle celle
            U(3:2:2*(N+1),1)=u0(x);
            
            % Creazione delle celle ghost con condizioni al bordo
            % periodiche
            U(1:2,1)=U(end-3:end-2,1);  % Cella ghost di sinistra
            U(end-1:end,1)=U(3:4,1);    % Cella ghost di destra
            

        case 3
            v0=@(x) 1.*(x<=0);
            w0=@(x) (0.5).*(x<=0) + 1.*(x>0);
            
            % Definizione delle approssimazioni iniziali date dal profilo
            % iniziale valutato nei centri delle celle
            U(3:2:2*(N+1),1)=v0(x);
            U(4:2:2*(N+1),1)=w0(x);

            % Creazione delle celle ghost con condizioni al bordo
            % di Free-Flow
            U(1:2,1)=U(3:4,1);                   % Cella ghost di sinistra
            U(end-1:end,1)=U(end-3:end-2,1);     % Cella ghost di destra
    end

end


% Flusso di Rusanov
function F=Rusanov_flux(Up,Um,lambda,A,N)
 
    % Creo la matrice per la vettorizzazione che sarà diagonale a blocchi
    % con A sulla diagonale
    blocks = repmat({A}, 1, N);    % Crea una cella con N copie di A
    K = blkdiag(blocks{:});        % Espande le celle e costruisce la matrice diagonale a blocchi
    
    F=0.5*K*(Up+Um)-(lambda/2)*(Up-Um);
end