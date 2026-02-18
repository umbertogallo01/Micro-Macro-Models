%% Foglio 1: Opinion Dynamics (Hegselmann-Krause model)

clc
clear
close all

disp("*** HEGSELMANN-KRAUSE MODEL ***")
N=input('Scegliere il numero di agenti del modello: ');
T=30;
h=0.05; % step-size

% Scelta del parametro epsilon
epsilon=input('Scegliere il parametro di bounded confidence parameter ε: ');

% Scelta del tipo di interazioni: symmetric or non-symmetric
scelta = menu('Scegliere il tipo di interazioni:', 'Simmetriche', 'Non Simmetriche');

% Scelta della distribuzione dei dati iniziali
distribution = menu('Scegliere la distribuzione iniziale degli agenti:', ...
    'U(0,1) : Distribuzione Uniforme', 'N(0.5,0.3^4) : Distribuzione Gaussiana Standard', ...
    'P(0.3,0.15,1,3) : Distribuzione di Pearson');

switch distribution
    case 1
        X_0=rand(N,1);                             % Distribuzione uniforme in [0,1]  

    case 2
        mu=0.5;
        sigma=0.3^4;
        X_0=mu + sigma * randn(N,1);                % Distribuzione Gaussiana

    case 3
        mu=0.3;
        sigma=0.15;
        skew=1;
        kurt=3;
        X_0 = pearsrnd(mu,sigma,skew,kurt,N,1);     % Distribuzione di Pearson
end


% Chiamata della funzione per calcolare l'evoluzione nel tempo delle
% opinioni
[X,m,t,type,m1,m2,Var,Ni,i]=HK(N,X_0,T,h,epsilon,scelta);

titolo="Modello di Hegselmann-Krause " + type;


% Visualizzazione dei risultati (SOLUZIONE NUMERICA + m1 + Varianza)

figure
sgtitle(titolo); % titolo generale
set(gcf, 'WindowState', 'maximized');   % Imposta la figura a schermo intero
plot(t,X,'LineWidth', 1);               % plotta le righe delle matrice X rispetto ai tempi t
%title("Evolution of Opinions in time")
xlabel('Tempo')
ylabel('x_i(t)')
xlim([0 t(size(X,2))])
axis square


figure
titolo2="Modello di Hegselmann-Krause " + type;
sgtitle(titolo2);
set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero
subplot(1,3,1)
plot(t,m1,'LineWidth', 2);
title('Media delle opinioni')
xlabel('Tempo')
ylabel('m1(t)')
ylim([0 1]);
xlim([0 t(size(X,2))])
axis square


subplot(1,3,3)
plot(t,Var,'LineWidth', 2);
title('Varianza delle opinioni')
xlabel('Tempo')
ylabel('Var(t)')
xlim([0 t(size(X,2))])
axis square

if scelta==2
    subplot(1,3,2)
    plot(t,m,'LineWidth',2);
    title('Moleplicità geometrica dell''atovalore 1 della matrice A')
    xlabel('Tempo')
    ylabel('m(t)')
    xlim([0 t(size(X,2))])
    axis square
end


figure
plot(t,Ni,'LineWidth',2);
set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero
title("|Ni(t,eps)| per l'agente i= " + num2str(i));
xlabel('Tempo')
ylabel('|Ni(t,eps)|')
xlim([0 t(size(X,2))])
axis square


% ANIMAZIONE
waitforbuttonpress;
figure

% Creazione grafico iniziale
pa = plot(t(1), X(:,1),'LineWidth', 1);
set(gcf, 'WindowState', 'maximized');  % Imposta la figura a schermo intero
titolo3="Evoluzione temporale delle opinioni nel modello di Hegselmann Krause " + type;
sgtitle(titolo3);
xlabel('Tempo')
ylabel('x_i(t)')
xlim([0 t(size(X,2))])
axis square; 

% Gestione dell'asse per limitare i dati
a = gca;
a.XLim = [min(t), t(size(X,2))];
a.YLim = [min(X(:)), max(X(:))];  % Imposta i limiti degli assi


% Ciclo per l'animazione
for n = 1:3:size(X,2)-1
    for i=1:N
        pa(i).XData=[pa(i).XData t(n+1)];
        pa(i).YData=[pa(i).YData X(i,n+1)];
    end

    % Rende visibile l'aggiornamento
    drawnow;
    pause(0.02)
end




%% Funzioni

function [X,m,t,type,m1,m2,Var,Ni,i]=HK(N,X_0,T,h,epsilon,scelta)
   
    % Inizializzazione dell'approssimazione numerica
    X=X_0;
    
    t=0;             % tempo iniziale
    tol=1.e-8;       % tolleranza


    % Inizializzazione di momento primo/secondo e Varianza
    m1=(1/N)*sum(X(:,1));
    m2=(1/N)*sum(X(:,1).^2);
    Var=m2-m1.^2;


    % Inizializzazione di un vettore che ha come componenti |Ni(t,eps)| :
    % consideriamo l'agente i=10
    i=10;
    Ni=[];


    % Calcolo della massima distanza tra gli agenti in modo da verificare
    % la condizioni sufficiente per avere un modello globale. Per farlo
    % costruiamo innanzitutto una matrice R che ha per colonne x_i ripetuto 
    % N volte e andiamo a sottrarre il vettore delle opinioni iniziali X0

    R = repmat(X_0', length(X_0), 1);
    D = abs(R - X_0);   
    max_dist = max(D(:));  % max su tutte le distanze - D(:) converte la matrice D in un vettore leggendo colonna per colonna    

    
    % Switch sulla scelta del tipo di interazioni
    switch scelta

        case 1 % Interazioni simmetriche
            
            m=-1; % Si intende che il valore di m non ci interessa

            if epsilon>max_dist % MODELLO GLOBALE
                
                type='Globale';

                % Preallochiamo la matrice di iterazione K
                A=(1/N)*ones(N);
                K=A-eye(N);
                

                n=1;
                
                % Al tempo 0, dato che siamo nel modello globale,
                % sicuramente Ni(0,eps)=N;
                Ni=N;

                while t(n)<T-h

                    % Applichiamo il metodo di Heun
                    X=[X X(:,n)+(h/2)*(K*X(:,n)+ K*(X(:,n)+h*K*X(:,n)))];

                    % Aggiorniamo il vettore dei tempi
                    t=[t t(n)+h];

                    n=n+1;

                    % Aggiorniamo il vettore Ni(t(n);eps))
                    Ni=[Ni N];

                    % Aggiornamento di momento primo/secondo e Varianza
                    m1=[m1 (1/N)*sum(X(:,n))];
                    m2=[m2 (1/N)*sum(X(:,n).^2)];
                    Var=[Var m2(n)-m1(n).^2];
                    

                    % CRITERIO D'ARRESTO
                    check=abs(Var(n-1)-Var(n));

                    if check < tol
                        return
                    end
                end
               

            else % MODELLO LOCALE SIMMETRICO

                type='Locale Simmetrico';
                n=1;

                while t(n)<T-h
                    
                    % Definizione di sigma_i
                    sigma=N;


                    % Costruzione della matrice K di iterazione
                    [~,I,K]=matrix_construction(X(:,n),epsilon,sigma);


                    % Aggiornamento del vettore Ni (ricorda i=10)
                    Ni=[Ni I(i,i)];

                    
                    % Applichiamo il metodo di Heun
                    X=[X X(:,n)+(h/2)*(K*X(:,n)+ K*(X(:,n)+h*K*X(:,n)))];

  
                    % Aggiorniamo il vettore dei tempi
                    t=[t t(n)+h];

                    n=n+1; 


                    % Aggiornamento di momento primo/secondo e Varianza
                    m1=[m1 (1/N)*sum(X(:,n))];
                    m2=[m2 (1/N)*sum(X(:,n).^2)];
                    Var=[Var m2(n)-m1(n).^2];
                        

                    % CRITERIO D'ARRESTO
                    check=abs(Var(n-1)-Var(n));

                    if check < tol
                        % Manca l'ultima iterazione per aggiornare il vettore Ni
                        [~,I,~]=matrix_construction(X(:,n),epsilon,sigma);
                        Ni=[Ni I(i,i)];
                        return
                    end
                end

                % Manca l'ultima iterazione per aggiornare il vettore Ni
                [~,I,~]=matrix_construction(X(:,n),epsilon,sigma);
                Ni=[Ni I(i,i)];

            end
        


        case 2 % INTERAZIONI NON SIMMETRICHE
            
            if epsilon>max_dist % Modello Globale
                error("C'e' un errore nella scelta di epsilon, perché se si considerano interazioni non simmetriche, " + ...
                    "allora epsilon non può essere più grande della massima distanza tra le opinioni al tempo 0");
            end

            type='Non simmetrico';
            n=1;
            m=[];

            while t(n)<T-h
                
                [A,K,sigma]=no_symmetric_model(X(:,n),epsilon,N);

                % Aggiornamento del vettore Ni (ricorda i=10)
                Ni=[Ni sigma(i)];


                % Salviamo la molteplicità geometrica dell'autovalore 1 di
                % A in un vettore m: per farlo usiamo il comando null(M)
                % che calcola una base ortonormale per il kernel della 
                % matrice M e successivamente usiamo il comando size(.,2)
                % che resituisce il numero di colonne della matrice 

                multiplicity=size(null(A-1*eye(size(A))),2);
                m =[m multiplicity];

                
               % Applichiamo il metodo di Heun
                X=[X X(:,n)+(h/2)*(K*X(:,n)+ K*(X(:,n)+h*K*X(:,n)))];


                % Aggiorniamo il vettore dei tempi
                t=[t t(n)+h];

                n=n+1; 
                

                % Aggiornamento di momento primo/secondo e Varianza
                m1=[m1 (1/N)*sum(X(:,n))];
                m2=[m2 (1/N)*sum(X(:,n).^2)];
                Var=[Var m2(n)-m1(n).^2];


                % CRITERIO D'ARRESTO
                check=abs(Var(n-1)-Var(n));
               
                if check < tol
                    % Manca l'ultima iterazione per aggiornare il vettore
                    % Ni e il vettore m
                    [A,~,sigma]=no_symmetric_model(X(:,n),epsilon,N);

                    Ni=[Ni sigma(i)];

                    multiplicity=size(null(A-1*eye(size(A))),2);
                    m =[m multiplicity];
                    return
                end            
            end

            % Manca l'ultima iterazione per aggiornare il vettore
            % Ni e il vettore m
            [A,~,sigma]=no_symmetric_model(X(:,n),epsilon,N);

            Ni=[Ni sigma(i)];

            multiplicity=size(null(A-1*eye(size(A))),2);
            m =[m multiplicity];
    end
end



% Funzione per la costruzione delle matrici A,I,K del modello locale
% simmetrico
function [A,I,K]=matrix_construction(X,epsilon,sigma)

    % 1) Calcolo della matrice delle distanze
    R = repmat(X', length(X), 1);
    D = abs(R - X);
    
    % 2) Calcolo una matrice NxN di 0/1 che individua tutti
    % gli elementi della matrice delle distanze D che sono
    % <=eps
    RI=(D<=epsilon);

    % 3) Calcolo della matrice A
    A=RI.*(1./sigma);

    % 4) Calcolo della matrice I~: sommiamo per righe gli
    % elementi della matrice RI booleana, in modo da
    % ottenere un vettore colonna che ha per componenti
    % |Ni(t,eps)|. Infine costruiamo una matrice che ha
    % questo vettore come diagonale
    I=diag(sum(RI,2));

    % 5) Calcolo della matrice K
    K=A-(1/sigma)*I;

end



% Funzione per la costruzione delle matrici A e K del modello locale
% NON-simmetrico e del vettore sigma
function [A,K,sigma]=no_symmetric_model(X,epsilon,N)

    % 1) Calcolo della matrice delle distanze
    R = repmat(X', length(X), 1);
    D = abs(R - X);
    
    % 2) Calcolo una matrice NxN di 0/1 che individua tutti
    % gli elementi della matrice delle distanze D che sono
    % <=eps
    RI=(D<=epsilon);

    % 3) Calcolo della matrice A
    sigma=sum(RI,2);
    A=RI.*(1./sigma);

    % 4) Calcolo della matrice K
    K=A-eye(N);

end


