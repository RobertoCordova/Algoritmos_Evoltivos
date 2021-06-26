%% Problema del agente viajero%%
clear;
close;
clc;
format long g;

figure
axis([0 1000 20000 45000])
xlabel('Generacion')
ylabel('Individuo mas apto')
title('Azul -TorneoK/L. Verde -Torneo P/U. Rojo -TorneoE P/U ')

format longG
%Cuerpo principal del programa%

%Cargar la tabla contenedora de distancios
load('tablaconfusion');
tablaconfusion = round(tablaconfusion);

%Creacion de la poblacion inicial
poblacion = CreacionPoblacion(500,18);

%Obtener aptitud
poblacion = Aptitud(poblacion,tablaconfusion);

%Ordenar por aptitud
poblacion = Ordenamiento(poblacion);

%Cruze PMx
nuevagen = pMx(poblacion,7,11);
nuevagen1 = pTm(poblacion,7,11);
nuevagen2 = pTe(poblacion,7,11,tablaconfusion);

%Obtener aptitud
nuevagen = Aptitud(nuevagen,tablaconfusion);
nuevagen1 = Aptitud(nuevagen1,tablaconfusion);
nuevagen2 = Aptitud(nuevagen2,tablaconfusion);

nuevagen = Ordenamiento(nuevagen);
nuevagen1 = Ordenamiento(nuevagen1);
nuevagen2 = Ordenamiento(nuevagen2);

for t = 1:1000
    hold on
    plot(t,nuevagen(1,19),'+blue',t,nuevagen1(1,19),'*green',...
         t,nuevagen2(1,19),'ored')
    pause(.01)
   
    %Cruze PMx
    nuevagen = pMx(nuevagen,7,11);
    nuevagen1 = pTm(nuevagen1,7,11);
    nuevagen2 = pTe(nuevagen2,7,11,tablaconfusion);
    
    nuevagen = Aptitud(nuevagen,tablaconfusion);
    nuevagen1 = Aptitud(nuevagen1,tablaconfusion);
    nuevagen2 = Aptitud(nuevagen2,tablaconfusion);
    
    nuevagen = Ordenamiento(nuevagen);
    nuevagen1 = Ordenamiento(nuevagen1);
    nuevagen2 = Ordenamiento(nuevagen2);
    
    %Mutacion
    nuevagen = MutacionN(nuevagen(:,1:size(nuevagen,2)-1));
    nuevagen1 = MutacionN(nuevagen1(:,1:size(nuevagen1,2)-1));
    nuevagen2 = MutacionN(nuevagen2(:,1:size(nuevagen2,2)-1));
    
    %Obtener aptitud
    nuevagen = Aptitud(nuevagen,tablaconfusion);
    nuevagen1 = Aptitud(nuevagen1,tablaconfusion);
    nuevagen2 = Aptitud(nuevagen2,tablaconfusion);
    
    %Ordenar por aptitud
    nuevagen = Ordenamiento(nuevagen);
    nuevagen1 = Ordenamiento(nuevagen1);
    nuevagen2 = Ordenamiento(nuevagen2);
end

%% Creacion de la poblacion %%
function pobl = CreacionPoblacion(ta,ciu)
in = 1:1:ciu;
pobl = zeros(ta,ciu);
for f = 1:ta
    for c = ciu:-1:1
        cambio = randi(c);
        guardar = in(1,c);
        in(1,c) = in(1,cambio);
        in(1,cambio) = guardar;
    end
    pobl(f,:) = in;
end
end

%% Aptitud de los individuos %%
function pobl = Aptitud(poblacion,tablaconfusion)
[fi,co] = size(poblacion);
distancia = zeros(fi,1);
for f = 1:fi
    for c = 1:co-1
        distancia(f,1) = distancia(f,1) + tablaconfusion(...
                                poblacion(f,c),poblacion(f,c+1));
    end
    distancia(f,1) = distancia(f,1) + tablaconfusion(...
                            poblacion(f,1),poblacion(f,c+1));
end
pobl = [poblacion(:,:) distancia(:,:)];
end

%% Ordenamiento %%
function pobl = Ordenamiento(poblacion)
[fi,co] = size(poblacion);
[in,k] = sort(poblacion(:,co),'ascend');
[fi,co] = size(poblacion);
pobl = zeros(fi,co);
for f = 1:fi
    pobl(f,:) = poblacion(k(f),:); 
end
end
function pobl = Ordenamiento1(poblacion)
[fi,co] = size(poblacion);
[in,k] = sort(poblacion(:,co),'descen');
[fi,co] = size(poblacion);
pobl = zeros(fi,co);
for f = 1:fi
    pobl(f,:) = poblacion(k(f),:); 
end
end

%% Seleccion torneo K/L %%
function [cont1,cont2] = TorneoKL(poblacion,contendientes)
f = 1;
[fi,co] = size(poblacion);
prosp1 = zeros(contendientes,co);
candidato = zeros(contendientes,1);
w = 0;
for t = 1:2
    while f <= contendientes
        candidato(f,1) = randi(fi);
        for che = 1:contendientes
            if candidato(f,1) ~= candidato(che,1)
                prosp1(f,:) = poblacion(candidato(f,1),:);
            else
                w = w + 1;
            end
        end
        if w <= 1
            f = f + 1;
        end
        w = 0;
    end
    prosp1 = Ordenamiento(prosp1);
    if t == 1
        cont1 = prosp1(1,1:co-1);
        f = 1;
    else
        cont2 = prosp1(1,1:co-1);
    end
end
end

%% Cruza por torneo K/L %%%
function nuevaGeneracion = pMx(poblacion,c1,c2)
[fi,co] = size(poblacion);
nuevaGeneracion = zeros(fi,co-1);
for p = 2:2:fi
    %Seleccion
    [padre,madre] = TorneoKL(poblacion,4);
    relac = zeros(2,c2-c1);
    relac(1,:) = padre(1,c1+1:c2);
    relac(2,:) = madre(1,c1+1:c2);
    padre(1,c1+1:c2) = relac(2,:);
    madre(1,c1+1:c2) = relac(1,:);
    for rr = 1:c2-c1
        for ch = 1:c2-c1
            for c = 1: length(padre)
                if (relac(2,ch) == padre(1,c)) && (c <= c1 || c > c2)
                    padre(1,c) = relac(1,ch);
                end
                if (relac(1,ch) == madre(1,c)) && (c <= c1 || c > c2)
                    madre(1,c) = relac(2,ch);
                end
            end
        end
    end
    checador(padre);
    nuevaGeneracion(p-1,:) = padre(1,:);
    nuevaGeneracion(p,:) = madre(1,:);
end
end

%% Mutacion %%
function poblacionn = MutacionN(poblacion)
poblacionn = poblacion;
p = 1; e = zeros(1,length(poblacionn)*.05);
while p <= size(e,2)
    t = 0;
    pe = randi([length(poblacion)-length(poblacionn)/4,length(poblacion)]);
    for b = 1:length(e)
        if pe == e(b)
            t = 1;
            break;
        end
    end
    if t == 0
        e(p) = pe;
        p = p + 1;
    end
end

for rr = 1:length(e)
    rango = randi([6,9]);
    zonac = randi(size(poblacionn,2)-rango);
    zonai = randi(size(poblacionn,2)-rango);
    el = poblacionn(e(rr),:);
    if (zonac+rango) < length(el)
        ety = [el(1:zonac), el((zonac+rango+1):length(el))];
    else
        ety =el(1:zonac);
    end
    extr = el((zonac+1):(zonac+rango));
    for cn = 1:length(extr)
        xc = randi(length(extr));
        axtr = extr(xc);
        while 1
            xz = randi(length(extr));
            if xc ~= xz
                ixtr = extr(xz);
                extr(xz) = axtr;
                extr(xc) = ixtr;
                break;
            end
        end
    end
    if (zonai+rango) < length(el)
        el = [ety(1:zonai),extr,ety((zonai+1):length(ety))];
        %el = [el(1:zonai),extr,el((zonai+rango+1):size(el,2))];
    else
        %error
        el = [ety(1:zonai),extr];
        %el = [el(1:zonai),extr];
    end
    poblacionn(e(rr),:) = el(:);
end
end

%% Cruza por torneo %%
function nuevag = pTm(p,c1,c2)
nuevag = zeros(size(p,1),size(p,2)-1);
for f = 0:((size(p,1)/2)-1)
    padre = p(f+1,1:size(p,2)-1);
    madre = p((size(p,1)-f),1:size(p,2)-1);
    relac = zeros(2,c2-c1);
    relac(1,:) = padre(1,c1+1:c2);
    relac(2,:) = madre(1,c1+1:c2);
    padre(1,c1+1:c2) = relac(2,:);
    madre(1,c1+1:c2) = relac(1,:);
    for rr = 1:c2-c1
        for ch = 1:c2-c1
            for c = 1: length(padre)
                if (relac(2,ch) == padre(1,c)) && (c <= c1 || c > c2)
                    padre(1,c) = relac(1,ch);
                end
                if (relac(1,ch) == madre(1,c)) && (c <= c1 || c > c2)
                    madre(1,c) = relac(2,ch);
                end
            end
        end
    end
    nuevag(f*2+1,:) = padre(1,:);
    nuevag(f*2+2,:) = madre(1,:);
end
end

%% Cruza por torneo con elitismo %%
function nuevag = pTe(p,c1,c2,tablaconfusion)
nuevag = zeros(size(p,1),size(p,2)-1);
for f = 0:((size(p,1)/2)-1)
    padreo = p(f+1,1:size(p,2)-1);
    madreo = p((size(p,1)-f),1:size(p,2)-1);    
    padre = p(f+1,1:size(p,2)-1);
    madre = p((size(p,1)-f),1:size(p,2)-1);
    relac = zeros(2,c2-c1);
    relac(1,:) = padre(1,c1+1:c2);
    relac(2,:) = madre(1,c1+1:c2);
    padre(1,c1+1:c2) = relac(2,:);
    madre(1,c1+1:c2) = relac(1,:);
    for rr = 1:c2-c1
        for ch = 1:c2-c1
            for c = 1: length(padre)
                if (relac(2,ch) == padre(1,c)) && (c <= c1 || c > c2)
                    padre(1,c) = relac(1,ch);
                end
                if (relac(1,ch) == madre(1,c)) && (c <= c1 || c > c2)
                    madre(1,c) = relac(2,ch);
                end
            end
        end
    end
    padreo = Aptitud(padreo,tablaconfusion);
    madreo = Aptitud(madreo,tablaconfusion);
    padre = Aptitud(padre,tablaconfusion);
    madre = Aptitud(madre,tablaconfusion);
    losmejore = Ordenamiento([padre;madre;padreo;madreo]);
    nuevag(f*2+1,:) = losmejore(1,1:size(losmejore,2)-1);
    nuevag(f*2+2,:) = losmejore(2,1:size(losmejore,2)-1);
end
end

%% Reviza que no exista repeticion %%
function checador(poblacion)
t = 0;
for n = 1:18
    for f = 1:size(poblacion,1)
        for c = 1:size(poblacion,2)
            if poblacion(f,c) == n
                t = t + 1;
            end
        end
        if t > 1
            disp("Error");
        end
        t = 0;
    end
end
end