%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Procesamiento de señales II %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Trabajo práctico final: Aplicación del Filtro de Kalman a %%%%%%%%
%%%%%%%%%%%%%% predicción de cotizaciones de commodities %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Autor: Gonzalo Quintana %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

SerieElegida=4;
%1) Cotización del Dólar en coronas suecas
%2) Cotización de la onza de oro en dólares
%3) Cotización de la tonelada de soja en dólares
%4) Cotización del Bitcoin en dólares
%5) Cotización del barril de petróleo en dólares
%6) Cotización del Dólar en pesos argentinos

DiasSinAct=15;
AnoSinAct=2015;
NPrediccionesExtra=30;
Nmontecarlo=1;
VectorN=[6, 4, 6, 4, 6, 10];%el 6 para el 1 ser funciona bien


switch SerieElegida
    case 1
        SerieTemporal=dlmread('Datos históricos USD SEK.txt');
    case 2
        SerieTemporal=dlmread('Datos históricos Futuros oro.txt');
    case 3
        SerieTemporal=dlmread('Datos históricos Futuros soja EE.UU.txt');
    case 4
        SerieTemporal=dlmread('Datos históricos BTC USD.txt');
        VectorQ=0.000001:0.000001:0.00005;
        VectorR=1:1:130;
        RvaracionqQ=50;
        QvaracionqR=0.00001;
    case 5
        SerieTemporal=dlmread('Datos históricos Futuros petróleo crudo WTI.txt');
        VectorQ=0.000001:0.00001:0.0005;
        VectorR=1:50;
        RvaracionqQ=10;
        QvaracionqR=0.0001;
    case 6
        SerieTemporal=dlmread('dolar.txt');
end

Nsistema=VectorN(SerieElegida);
Ncoeficientes=Nsistema*(Nsistema+1)/2;
N=Nsistema+Ncoeficientes;

if SerieElegida~=6
    SerieTemporal=flip(SerieTemporal);
end

ECM=zeros(size(VectorQ));
ECMLineal=zeros(size(VectorQ));

for iteracionQ=1:length(VectorQ)
    
    
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
    QkLineal=VectorQ(iteracionQ);
    RkLineal=RvaracionqQ;

    XpredichoLineal=zeros(2,SerieTemporal(end,1)-SerieTemporal(1,1)+1);
    XLineal=zeros(2,SerieTemporal(end,1)-SerieTemporal(1,1)+1);
    TiempoDiscretoLineal=TiempoDiscreto;
    ErrorLineal=zeros(1,length(SerieTemporal));
    TiemposSinActualizacion=(AnoSinAct-1900)*365:(AnoSinAct-1900)*365+DiasSinAct;
    TiemposSinActualizacion=[TiemposSinActualizacion, TiemposSinActualizacion(end)];

    XLineal(:,1)=zeros(2,1);
    PkLineal=eye(2,2);

    %Matrices del sistema al usar el modelo no lineal
    A11=zeros(Nsistema,Nsistema);
    A12=zeros(Nsistema,Ncoeficientes);
    A21=zeros(Ncoeficientes,Nsistema);
    A22=eye(Ncoeficientes,Ncoeficientes);

    Bk=ones(N,1);
    Ck=[1,zeros(1,N-1)];
    Qk=VectorQ(iteracionQ); %covarianza del ruido de proceso

    Rk=RvaracionqQ; %covarianza del ruido de medición

    %Condición inicial para el filtro
    X(1:Nsistema,1)=ones(Nsistema,1);
    X(Nsistema+1:end,1)=zeros(Ncoeficientes,1);
    Pk=eye(N,N);

    j=2;
    tsa=1;
    for k=2:TiemposSinActualizacion(end)+10-SerieTemporal(1,1)

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

    %predicción modelo no lineal
        Xpredicho(:,k)=Ak*X(:,k-1);
        Ppredicho=Ak*Pk*Ak'+Bk*Qk*Bk';

    %predicción modelo lineal
        XpredichoLineal(:,k)=AkLineal*XLineal(:,k-1);
        PpredichoLineal=AkLineal*PkLineal*AkLineal'+BkLineal*QkLineal*BkLineal';



        if TiempoDiscreto(k)==TiemposSinActualizacion(tsa)

                        %modelo no lineal
                        X(:,k)=Xpredicho(:,k);
                        Pk=Ppredicho;

                        %modelo lineal
                        XLineal(:,k)=XpredichoLineal(:,k);
                        PkLineal=PpredichoLineal;

                        if TiempoDiscreto(k)==SerieTemporal(j,1)
                            Error(j)=abs(X(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);
                            ErrorLineal(j)=abs(XLineal(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);
                        end

                        if tsa==1
                            jInicial=j;
                        end
                        if tsa==length(TiemposSinActualizacion)-1
                            jFinal=j;
                        end

                        tsa=tsa+1;
                        j=j+1;
        else

                if TiempoDiscreto(k)==SerieTemporal(j,1)

                    yk=SerieTemporal(j,2);

                    %actualización modelo lineal
                    KkLineal=PpredichoLineal*CkLineal'/(RkLineal+CkLineal*PpredichoLineal*CkLineal');
                    XLineal(:,k)=XpredichoLineal(:,k)+KkLineal*(yk-CkLineal*XpredichoLineal(:,k));
                    PkLineal=(eye(2,2)-KkLineal*CkLineal)*PpredichoLineal;
                    ErrorLineal(j)=abs(XLineal(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);

                    %actualización modelo no lineal
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
    end

    ECM(iteracionQ)=sum(Error(jInicial:jFinal))/length(TiemposSinActualizacion);
    ECMLineal(iteracionQ)=sum(ErrorLineal(jInicial:jFinal))/length(TiemposSinActualizacion);

end



TiempoReescalado=TiempoDiscreto/365+1900;
Anios=round(TiempoReescalado(1)):round(TiempoReescalado(end));

figure(1)
plot(VectorQ,ECM);
hold on;
plot(VectorQ,ECMLineal);
legend('ECM','ECM Lineal','Interpreter','Latex');
grid on;
xlabel('Varianza del ruido de proceso','Interpreter','Latex');
ylabel('ECM','Interpreter','Latex');
%ylim([0 0.2]);
print('ECMvsQ_2015.png','-dpng');

for iteracionR=1:length(VectorR)
    
    
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
    QkLineal=QvaracionqR;
    RkLineal=VectorR(iteracionR);

    XpredichoLineal=zeros(2,SerieTemporal(end,1)-SerieTemporal(1,1)+1);
    XLineal=zeros(2,SerieTemporal(end,1)-SerieTemporal(1,1)+1);
    TiempoDiscretoLineal=TiempoDiscreto;
    ErrorLineal=zeros(1,length(SerieTemporal));
    TiemposSinActualizacion=(AnoSinAct-1900)*365:(AnoSinAct-1900)*365+DiasSinAct;
    TiemposSinActualizacion=[TiemposSinActualizacion, TiemposSinActualizacion(end)];

    XLineal(:,1)=zeros(2,1);
    PkLineal=eye(2,2);

    %Matrices del sistema al usar el modelo no lineal
    A11=zeros(Nsistema,Nsistema);
    A12=zeros(Nsistema,Ncoeficientes);
    A21=zeros(Ncoeficientes,Nsistema);
    A22=eye(Ncoeficientes,Ncoeficientes);

    Bk=ones(N,1);
    Ck=[1,zeros(1,N-1)];
    Qk=QvaracionqR; %covarianza del ruido de proceso

    Rk=VectorR(iteracionR); %covarianza del ruido de medición

    %Condición inicial para el filtro
    X(1:Nsistema,1)=ones(Nsistema,1);
    X(Nsistema+1:end,1)=zeros(Ncoeficientes,1);
    Pk=eye(N,N);

    j=2;
    tsa=1;
    for k=2:TiemposSinActualizacion(end)+10-SerieTemporal(1,1)

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

    %predicción modelo no lineal
        Xpredicho(:,k)=Ak*X(:,k-1);
        Ppredicho=Ak*Pk*Ak'+Bk*Qk*Bk';

    %predicción modelo lineal
        XpredichoLineal(:,k)=AkLineal*XLineal(:,k-1);
        PpredichoLineal=AkLineal*PkLineal*AkLineal'+BkLineal*QkLineal*BkLineal';



        if TiempoDiscreto(k)==TiemposSinActualizacion(tsa)

                        %modelo no lineal
                        X(:,k)=Xpredicho(:,k);
                        Pk=Ppredicho;

                        %modelo lineal
                        XLineal(:,k)=XpredichoLineal(:,k);
                        PkLineal=PpredichoLineal;

                        if TiempoDiscreto(k)==SerieTemporal(j,1)
                            Error(j)=abs(X(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);
                            ErrorLineal(j)=abs(XLineal(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);
                        end

                        if tsa==1
                            jInicial=j;
                        end
                        if tsa==length(TiemposSinActualizacion)-1
                            jFinal=j;
                        end

                        tsa=tsa+1;
                        j=j+1;
        else

                if TiempoDiscreto(k)==SerieTemporal(j,1)

                    yk=SerieTemporal(j,2);

                    %actualización modelo lineal
                    KkLineal=PpredichoLineal*CkLineal'/(RkLineal+CkLineal*PpredichoLineal*CkLineal');
                    XLineal(:,k)=XpredichoLineal(:,k)+KkLineal*(yk-CkLineal*XpredichoLineal(:,k));
                    PkLineal=(eye(2,2)-KkLineal*CkLineal)*PpredichoLineal;
                    ErrorLineal(j)=abs(XLineal(1,k)-SerieTemporal(j,2))./SerieTemporal(j,2);

                    %actualización modelo no lineal
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
    end

    ECM(iteracionR)=sum(Error(jInicial:jFinal))/length(TiemposSinActualizacion);
    ECMLineal(iteracionR)=sum(ErrorLineal(jInicial:jFinal))/length(TiemposSinActualizacion);

end

figure(2)
plot(VectorR,ECM);
hold on;
plot(VectorR,ECMLineal);
legend('ECM','ECM Lineal','Interpreter','Latex');
grid on;
xlabel('Varianza del ruido de proceso','Interpreter','Latex');
ylabel('ECM','Interpreter','Latex');
xlim([0 130]);
print('ECMvsR_2015.png','-dpng');