%% ES 2 Foglio 3- Linearized Euler Equations
% Write a code which numerically solves (5), endowed with an initial 
% condition U(x, t = 0) = U0(x), with an explicit first order finite volume
% scheme. In particular, consider the explicit Euler method for the time 
% discretization combined with the Rusanov numerical flux computed on 
% piecewise constant reconstruction of the solution.

clc
clear
close all

disp("*** Linearized Euler Equations ***");

% Parametri fissati per lo schema numerico
pho=1;
v=0;
gamma=1.4;
p=1;
A=[v pho 0; 0 v 1/pho ; 0 gamma*p v];
xm=-1;
xM=1;
T=0.5;
N=200;                  % numero di celle
dx=(xM-xm)/N;
x=xm+((1:N)-0.5).*dx;   % centri delle celle


% Calcolo del time-step con le condizioni di stabilità CFL
lambda=max(abs(eig(A)));
C=0.5; % costante di sicurezza
dt=(dx/lambda) * C;
t = 0:dt:T;
if t(end) < T           % non ricopro [0,Tfin] con dt della CFL 
    Nt = length(t)+1;

else                    % ricoprio [0,Tfin] con dt della CFL
    Nt = length(t);
end
dt = T/(Nt-1);
t=0:dt:T;


% Inizializzazione dell'approssimazione numerica
U=initial_conditions(x,N,Nt);


% Schema a volume finito
n=1;

while t(n)<T

    % Calcolo dei flussi di Rusanov
    Fd=Rusanov_flux(U(7:3*(N+2),n),U(4:3*(N+1),n),lambda,A,N);
    Fs=Rusanov_flux(U(4:3*(N+1),n),U(1:3*N,n),lambda,A,N);
    
    % Schema di Eulero esplicito
    U(4:3*(N+1),n+1) = U(4:3*(N+1),n) - (dt/dx)*(Fd-Fs);
   
    % Riapplichiamo le condizioni al bordo di free-flow
    U(1:3,n+1)=U(4:6,n+1);                  % cella ghost di sinistra
    U(end-2:end,n+1)=U(end-5:end-3,n+1);    % cella ghost di destra
      
    n=n+1;
end


% Visualizzazione grafica dei risultati

% Animazione
f1=plot(x,U(4:3:3*(N+1),1),'r',x,U(5:3:3*(N+1),1),'g',x,U(6:3:3*(N+1),1),'b');
set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero
grid on
axis square
xlabel('X')
ylabel('U(x,t)')
legend('ρ(x,t)','v(x,t)','p(x,t)')
title(['Tempo = ', num2str(0)]);
ylim([min(U,[],"all")-1 max(U,[],"all")+1])

waitforbuttonpress
% Ciclo per l'animazione
for i = 1:n-1
    f1(1).YData=U(4:3:3*(N+1),i+1);
    f1(2).YData=U(5:3:3*(N+1),i+1);
    f1(3).YData=U(6:3:3*(N+1),i+1);

    % Rende visibile l'aggiornamento
    drawnow;

    % aggiorna solo il titolo dell’asse della figura animata
    title(['Tempo = ', num2str(t(i+1))]);
    pause(0.001);
end


figure(2)
i=find(t==0.5); % Individuo l'indice del vettore dei tempi corrispondente a t=0.5
plot(x,U(4:3:3*(N+1),i),'r',x,U(5:3:3*(N+1),i),'g',x,U(6:3:3*(N+1),i),'b');
set(gcf, 'WindowState', 'maximized');  
axis square
grid on
xlabel('x')
ylabel('U(x,t=0.5)')
legend('ρ(x,t)','v(x,t)','p(x,t)')
title(['Soluzione numerica per t= ', num2str(0.5)]);
ylim([min(U,[],"all")-1 max(U,[],"all")+1])


%% Funzioni

% Condizioni iniziali
function U=initial_conditions(x,N,Nt)

    %Inizializzazione dell'approssimazione numerica con le due celle ghost
    U=zeros(3*(N+2),Nt);
    
    rho0=@(x) 1.*(x<=0) + (0.2).*(x>0);
    v0=@(x) 0;
    p0=@(x) 1.*(x<=0) + (0.2).*(x>0);

    % Definizione delle approssimazioni iniziali date dal profilo
    % iniziale valutato nei centri delle celle
    U(4:3:3*(N+1),1)=rho0(x);
    U(5:3:3*(N+1),1)=v0(x);
    U(6:3:3*(N+1),1)=p0(x);

    % Creazione delle celle ghost con condizioni al bordo
    % di Free-Flow
    U(1:3,1)=U(4:6,1);                  % cella ghost di sinistra
    U(end-2:end,1)=U(end-5:end-3,1);    % cella ghost di destra
end



% Flusso di Rusanov
function F=Rusanov_flux(Up,Um,lambda,A,N)
 
    % Creo la matrice per la vettorizzazione che sarà diagonale a blocchi
    % con A sulla diagonale
    blocks = repmat({A}, 1, N);  % Crea una cella con N copie di A
    K = blkdiag(blocks{:});      % Espande le celle e costruisce la matrice a blocchi
    
    F=0.5*K*(Up+Um)-(lambda/2)*(Up-Um);
end