%   Este script obtiene los coeficientes de la serie de Fourier, su
%   espectro y su representacion grafica a partir de vectores de
%   tiempo y datos de una onda periodica.
%
%   FFT_Series_Fourier_Coeficients.m
%
%   2015-11-16 por Raul Hurtado Gavilan

close all
clc
clear
format long

%Lectura de los datos y vectores de datos
Datos = csvread('Datos_ondas_ensayo.csv');
frecuencia = 50; %Frecuencia original de la onda en hertz
w = 2*pi*frecuencia; %Frecuencia natural en rad/s
t = Datos(:,5); t=t(1:length(t)/2); %Vector de tiempo de la funcion original
f = Datos(:,6); f=f(1:length(f)/2); %Vector de ordenadas de la funcion original

N = length(t); %Numero de elementos de los vectores de datos

%Calculo de la transformada rapida de Fourier para obtener los coeficientes
%de las series de Fourier de la funcion reconstruida
F = fft(f,N); %Transformada rapida de Fourier

F = F(1:N/2); %La transformada rapida de Fourier espeja los resultados, 
              %por tanto solo son utiles la primera mitad
n = length(F);

A = 2*real(F(2:n))/(N); %Vector de coeficientes de los cosenos (de 2 a n) 
                        %de la series de Fourier
A0= 2*real(F(1))/N; %Elemento a0 de la series de Fourier
B =-2*imag(F(2:n))/(N); %Vector de coeficientes de los senos (de 2 a n) 
                         %de la series de Fourier

t_rec=t(1:n-1); %Vector de tiempo con espaciado doble al original
                %para poder representar la funcion reconstruida
a=0; %Contador auxiliar
for i=1:n-1 %Bucle que hace el nuevo vector de tiempo
    t_rec(i) = t(i+a);
    a = a + 1;
end

C=zeros(n-1,1); %Bucle para vector de sumatorio de los cosenos con sus coeficientes
for i=1:n-1
    for j=1:n-2
       C(i) = C(i) + A(j)*cos((j)*w*t_rec(i));
    end
end
 
S=zeros(n-1,1); %Bucle para vector de sumatorio de los senos con sus coeficientes
for i=1:n-1
    for j=1:n-2
        S(i) = S(i) + B(j)*sin((j)*w*t_rec(i));
    end
end
 
for i=1:n-1 %Bucle para vector de la funcion reconstruida
    Y(i) = A0 + C(i) + S(i);
end

%Grafica de la funcion original enfrentada a la funcion reconstruida
figure(1), plot(t(1:N-2),f(1:N-2),'r'), hold on, plot(t_rec,Y,'b'),
title('Funcion original vs funcion reconstruida'),
xlabel('Tiempo en segundos'), ylabel('Tension en voltios'),
legend('Funcion original','Funcion reconstruida')

%Funcion reconstruida y sus armonicos fundamentales hasta el septimo
figure(2), plot(t(1:N-2),f(1:N-2),'color',[1,0.6,0.2]), hold all,
for i=1:8
    plot(t_rec, A0 + A(i)*cos(i*w*t_rec) + B(i)*sin(i*w*t_rec))
end
title('Funcion original y sus primeros armonicos'),
xlabel('Tiempo en segundos'), ylabel('Tension en voltios'),
S_legend = legend('Funcion original','Armonico 1','Armonico 2','Armonico 3'...
    ,'Armonico 4','Armonico 5','Armonico 6','Armonico 7');
set(S_legend,'FontSize',5);

%Espectro de frecuencias
ms = frecuencia*N/2; %Muestras por segundo
phz = ms/n; %Paso en Hertz 
hertz = 0:phz:ms-phz; %Vector de frencuecia en Hertz
figure(3), plot(hertz(1:20),abs(F(1:20))/n) %Grafica de magnitud-frecuencia
title('Espectro de frecuencias'), xlabel('Hertz'), ylabel('Magnitud')
