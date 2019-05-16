%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Procesamiento de se�ales II %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Trabajo pr�ctico final: Aplicaci�n del Filtro de Kalman a %%%%%%%%
%%%%%%%%%%%%%% predicci�n de cotizaciones de commodities %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Autor: Gonzalo Quintana %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

SerieElegida=4;
%1) Cotizaci�n del D�lar en coronas suecas
%2) Cotizaci�n de la onza de oro en d�lares
%3) Cotizaci�n de la tonelada de soja en d�lares
%4) Cotizaci�n del Bitcoin en d�lares
%5) Cotizaci�n del barril de petr�leo en d�lares
%6) Cotizaci�n del D�lar en pesos argentinos

NPrediccionesExtra=0;
Nmontecarlo=1;
VectorN=[6, 4, 6, 4, 3, 10];%el 6 para el 1 ser funciona bien

VectorQk=[0.001, 0.0001, 0.00001, 0.00001, 0.0001, 0.0001];
VectorRk=[1, 100, 100, 50, 10, 1];

VectorQkLineal=[1, 1, 1, 1, 1, 1];
VectorRkLineal=[1, 1, 1, 1, 1, 1];

switch SerieElegida
    case 1
        SerieTemporal=dlmread('Datos hist�ricos USD SEK.txt');
    case 2
        SerieTemporal=dlmread('Datos hist�ricos Futuros oro.txt');
    case 3
        SerieTemporal=dlmread('Datos hist�ricos Futuros soja EE.UU.txt');
    case 4
        SerieTemporal=dlmread('Datos hist�ricos BTC USD.txt');
    case 5
        SerieTemporal=dlmread('Datos hist�ricos Futuros petr�leo crudo WTI.txt');
    case 6
        SerieTemporal=dlmread('dolar.txt');
end

Nsistema=VectorN(SerieElegida);
Ncoeficientes=Nsistema*(Nsistema+1)/2;
N=Nsistema+Ncoeficientes;

if SerieElegida~=6
    SerieTemporal=flip(SerieTemporal);
end
Xpredicho=zeros(N,SerieTemporal(end,1)-SerieTemporal(1,1)+1);
X=zeros(size(Xpredicho));
TiempoDiscreto=SerieTemporal(1,1):SerieTemporal(end,1);
LongSerie=length(Xpredicho);
Error=zeros(1,length(SerieTemporal));


%Matrices del sistema al usar el modelp lineal
T=1;
AkLineal=[1, T;
    0, 1];
BkLineal=[1/2*T^2;
    T;];
CkLineal=[1, 0];
QkLineal=VectorQkLineal(SerieElegida);
RkLineal=VectorRkLineal(SerieElegida);

XpredichoLineal=zeros(2,SerieTemporal(end,1)-SerieTemporal(1,1)+1);
XLineal=zeros(2,SerieTemporal(end,1)-SerieTemporal(1,1)+1);
TiempoDiscretoLineal=TiempoDiscreto;
ErrorLineal=zeros(1,length(SerieTemporal));

XLineal(:,1)=zeros(2,1);
PkLineal=eye(2,2);

%Matrices del sistema al usar el modelo no lineal
A11=zeros(Nsistema,Nsistema);
A12=zeros(Nsistema,Ncoeficientes);
A21=zeros(Ncoeficientes,Nsistema);
A22=eye(Ncoeficientes,Ncoeficientes);

Bk=ones(N,1);
Ck=[1,zeros(1,N-1)];
Qk=VectorQk(SerieElegida); %covarianza del ruido de proceso
          
Rk=VectorRk(SerieElegida); %covarianza del ruido de medici�n

%Condici�n inicial para el filtro
X(1:Nsistema,1)=ones(Nsistema,1);
X(Nsistema+1:end,1)=zeros(Ncoeficientes,1);
Pk=eye(N,N);

j=2;
for k=2:LongSerie
    
%Matriz Ak para el modelo no lineal
    cursor=Nsistema+1;
    for i=1:Nsistema
        ceros=zeros(1,i-1);
        A11(i,:)=[ceros , X(cursor:cursor+Nsistema-1-length(ceros),k-1)'];
        cursor=cursor+Nsistema-1-length(ceros)+1;
        ceros2=zeros(1,(Nsistema-i)*((Nsistema-i)+1)/2);
        ceros1=zeros(1,Ncoeficientes-length(ceros2)-(Nsistema-i+1));
        A12(i,:)=[ceros1,X(i:Nsistema,k-1)',ceros2];
    end     
    
    Ak=[A11,A12;
       A21,A22;];

%predicci�n modelo no lineal
    Xpredicho(:,k)=Ak*X(:,k-1);
    Ppredicho=Ak*Pk*Ak'+Bk*Qk*Bk';
    
%predicci�n modelo lineal
    XpredichoLineal(:,k)=AkLineal*XLineal(:,k-1);
    PpredichoLineal=AkLineal*PkLineal*AkLineal'+BkLineal*QkLineal*BkLineal';
    
    if TiempoDiscreto(k)==SerieTemporal(j,1)
        yk=SerieTemporal(j,2);
        
        %actualizaci�n modelo lineal
        KkLineal=PpredichoLineal*CkLineal'/(RkLineal+CkLineal*PpredichoLineal*CkLineal');
        XLineal(:,k)=XpredichoLineal(:,k)+KkLineal*(yk-CkLineal*XpredichoLineal(:,k));
        PkLineal=(eye(2,2)-KkLineal*CkLineal)*PpredichoLineal;
        ErrorLineal(j)=abs(XLineal(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);
        
        %actualizaci�n modelo no lineal
        Kk=Ppredicho*Ck'/(Rk+Ck*Ppredicho*Ck');
        X(:,k)=Xpredicho(:,k)+Kk*(yk-Ck*Xpredicho(:,k));
        Pk=(eye(N,N)-Kk*Ck)*Ppredicho;
        Error(j)=abs(X(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);
        j=j+1;
        
    else    
        %modelo no lineal
        X(:,k)=Xpredicho(:,k);
        Pk=Ppredicho;
        
        %modelo lineal
        XLineal(:,k)=XpredichoLineal(:,k);
        PkLineal=PpredichoLineal;
    end
end

Xpredicho=[Xpredicho, zeros(N,NPrediccionesExtra)];
X=[X, zeros(N,NPrediccionesExtra)];

XpredichoLineal=[XpredichoLineal, zeros(2,NPrediccionesExtra)];
XLineal=[XLineal, zeros(2,NPrediccionesExtra)];

%Finalmente, se realizan predicciones acerca del comportamiento de la serie
%temporal elegida en base a la informaci�n aprendida anteriormente
for iteraciones=1:Nmontecarlo
    
    Xaux=zeros(N,NPrediccionesExtra+1);
    XpredichoAux=zeros(N,NPrediccionesExtra+1);
    XlinealAux=zeros(2,NPrediccionesExtra+1);
    XpredichoLinealAux=zeros(2,NPrediccionesExtra+1);

    Xaux(:,1)=X(:,LongSerie);
    XpredichoAux(:,1)=Xpredicho(:,LongSerie);
    XLinealAux(:,1)=XLineal(:,LongSerie);
    XpredichoLinealAux(:,1)=XpredichoLineal(:,LongSerie);

    for k=2:NPrediccionesExtra+1

        cursor=Nsistema+1;
        for i=1:Nsistema
            ceros=zeros(1,i-1);
            A11(i,:)=[ceros , Xaux(cursor:cursor+Nsistema-1-length(ceros),k-1)'];
            cursor=cursor+Nsistema-1-length(ceros)+1;
            ceros2=zeros(1,(Nsistema-i)*((Nsistema-i)+1)/2);
            ceros1=zeros(1,Ncoeficientes-length(ceros2)-(Nsistema-i+1));
            A12(i,:)=[ceros1,Xaux(i:Nsistema,k-1)',ceros2];
        end     

        Ak=[A11,A12;
           A21,A22;];

    %predicci�n no lineal
        XpredichoAux(:,k)=Ak*Xaux(:,k-1);%+Qk*randn;
        Ppredicho=Ak*Pk*Ak'+Bk*Qk*Bk';

    %predicci�n lineal
        XpredichoLinealAux(:,k)=AkLineal*XLineal(:,k-1);%+QkLineal*randn;
        PpredichoLineal=AkLineal*PkLineal*AkLineal'+BkLineal*QkLineal*BkLineal';

    %No hay actualizaci�n
        Xaux(:,k)=XpredichoAux(:,k); 
        XLinealAux(:,k)=XpredichoLinealAux(:,k);

    end
    
    Xpredicho(:,LongSerie+1:end)=Xpredicho(:,LongSerie+1:end)+XpredichoAux(:,2:end);
    X(:,LongSerie+1:end)=X(:,LongSerie+1:end)+Xaux(:,2:end);
    XpredichoLineal(:,LongSerie+1:end)=XpredichoLineal(:,LongSerie+1:end)+XpredichoLinealAux(:,2:end);
    XLineal(:,LongSerie+1:end)=XLineal(:,LongSerie+1:end)+XLinealAux(:,2:end);
end

Xpredicho(:,LongSerie+1:end)=Xpredicho(:,LongSerie+1:end)/Nmontecarlo;
X(:,LongSerie+1:end)=X(:,LongSerie+1:end)/Nmontecarlo;
XpredichoLineal(:,LongSerie+1:end)=XpredichoLineal(:,LongSerie+1:end)/Nmontecarlo;
XLineal(:,LongSerie+1:end)=XLineal(:,LongSerie+1:end)/Nmontecarlo;

TiempoExtra=TiempoDiscreto(end)+1:TiempoDiscreto(end)+NPrediccionesExtra;
TiempoDiscreto=[TiempoDiscreto,TiempoExtra];
TiempoReescalado=TiempoDiscreto/365+1900;
Anios=round(TiempoReescalado(1)):round(TiempoReescalado(end));

figure(1)
plot(TiempoReescalado,X(1,:)); %Grafico X en lugar de las predicciones para un mejor an�lisis
hold on;
plot(TiempoReescalado,XLineal(1,:));
hold on;
plot(SerieTemporal(:,1)/365+1900,SerieTemporal(:,2));
legend('Valores predichos por el modelo no lineal','Valores predichos por el modelo lineal','Valores reales','Interpreter','Latex');

switch SerieElegida
    case 1
        title('Cotizaci�n del D�lar en Coronas Suecas');
    case 2
        title('Cotizaci�n de la onza de oro en d�lares');
    case 3
        title('Cotizaci�n de la tonelada de soja en d�lares');
    case 4
        title('Cotizaci�n del Bitcoin en d�lares');
    case 5
        title('Cotizaci�n del barril de petr�leo en d�lares');
    case 6
        title('Cotizaci�n del Dolar en pesos argentinos');
end

grid on;
MaxVal=1.1*max(SerieTemporal(:,2));
ylim([0 MaxVal]);
xlim([TiempoReescalado(1) TiempoReescalado(end)]);
set(gca,'xtick',Anios,'XTickLabelRotation',45);
%print('BitcoinSerie.png','-dpng');

%Se observa que la estimaci�n con el modelo alineal tiene menos picos 
%(con QkLineal=1). Al bajar el valor de QkLineal, la estimaci�n del modelo 
%lineal se hace incluso m�s suave que la del modelo alineal pero tiene un 
%retardo mayor(ya con Qk=0.01 se nota este efecto). Al aumentar QkLineal 
%por encima de la unidad, los picos son cada vez mayores, aunque el retardo
%es menor.

figure(2)
for i=Nsistema+1:N
    plot(TiempoReescalado,Xpredicho(i,:));
    hold on;
end
%ylim([0 1]);
xlim([TiempoReescalado(1) TiempoReescalado(end)]);
ylim([0 1]);
set(gca,'xtick',Anios,'XTickLabelRotation',45);
%title('Evoluci�n de los coeficientes del modelo AR');
grid on;
%print('BitcoinCoeficientes.png','-dpng');

figure(3)
plot(SerieTemporal(:,1)/365+1900,Error);
hold on;
plot(SerieTemporal(:,1)/365+1900,ErrorLineal);
grid on;
set(gca,'xtick',Anios,'XTickLabelRotation',45);
%title('Evoluci�n del error de estimaci�n');
legend('Modelo no lineal','Modelo lineal');
ylim([0 0.1]);
