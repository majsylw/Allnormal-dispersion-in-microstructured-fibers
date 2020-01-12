%Symulacja ukazuj¹ca proces tworzenia solitonu rzêdu 3.
%Dudley, J. M., G. Genty oraz S. Coen, 
%„Supercontinuum generation in photoniccrystal fiber,” 
%Review of Modern Physcis, 78(4), str. 1135–1184, pazdziernik 2006.
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.

clear all;
clc;
format long e

%--------------------------------------------------------------------------
%Parametry numeryczne
%--------------------------------------------------------------------------

T = 12.5;                   %szerokoœæ okna czasowego [ps]
N = 2^13;                   %iloœæ punktów siatki czasowej
dt = T/N;                   %krok czasowy [ps]
dw = 1/T;                   %krok czêstotliwoœciowy [THz]

t = zeros(1,N);             %siatka punktów czasowych
w = zeros(1,N);             %siatka punktów czêstotliwoœciowych
for i = 1 : N
    t(i) = -T/2 + i*dt;
    w(i) = 2 * pi * (-N/2 + i) * dw;
end

%--------------------------------------------------------------------------
% Warunki pocz¹tkowe - impuls pocz¹tkowy
%--------------------------------------------------------------------------

%Af = zeros(M,N);            %A(z,w)
%At = zeros(M,N);            %A(z,t)

P0 = 1250;                  %maksimum mocy impulsu pocz¹tkowego [W]
T0 = 0.0284;                %czas trwania impulsu [ps]

%Impuls pocz¹tkowy w domenie czasowej (sekans hiperboliczny)
At(1,:) = sqrt(P0) .* sech(t./T0);

%Impuls pocz¹tkowy w domenie czêstotliwoœciowej
Af(1,:) = FT_scal(t,At(1,:));


%--------------------------------------------------------------------------
%Parametry œwiat³owodu
%--------------------------------------------------------------------------
Beta = [0 -11.83];                           %pochodne beta(i) [ps^i/km]
m = length(Beta);                            %liczba pochodnych beta
alpha = 0;                                   %wspó³czynnik t³umiennoœci [W/km]
gamma = 110;                                 %wspó³czynnik nieliniowoœci [1/W/km]

fr = 0;                                      %
tau1 = 0.0122;                               %[ps]
tau2 = 0.032;                                %[ps]

c = 2.99792458*1e-7;                         %prêdkoœæ œwiat³a w pró¿ni [km/ps]
lambda_centralna = 835*1e-12;                %centralna d³ugoœæ fali [km]
w_c = (2.0*pi*c)/lambda_centralna;           %czêstoœæ centralna [2*pi*THz]

Ld = (T0^2)/(abs(Beta(2)));                  %droga dyspersji [km]
Ln = 1/(P0 * gamma);                         %droga nieliniowoœci [km]

z = 2 * (pi/2) * Ld;                         %d³ugoœæ œwiat³owodu [km]
dz = 5e-7;                                	 %pocz¹tkowy krok po³o¿eniowy [km]

%--------------------------------------------------------------------------
%Operator dyspersji
%--------------------------------------------------------------------------

B = 0;                      %Beta(w)
for k = 2:m
    B = B + Beta(k)/factorial(k) * w.^k;
end

D = 1i * B - alpha/2;       %operator liniowy

%--------------------------------------------------------------------------
%Rozpraszanie Ramana
%--------------------------------------------------------------------------

hR_t = (tau1^2+tau2^2)/tau1/tau2^2*exp(-t./tau2).*sin(t./tau1);
hR_t(t<0) = 0;
hR_f = fftshift(ifft(fftshift(hR_t))) * dt * N;

%--------------------------------------------------------------------------
%RK4IP
%--------------------------------------------------------------------------
ii = 2;
len = 0;                %siatka punktów po³o¿eniowych

propagedlength = 0;
fprintf(1, '\nStart...      ');
tic

while propagedlength < z,
    
    if (dz + propagedlength) > z
        dz = z - propagedlength;
    end
   
    AI = exp(D.*dz/2.0) .* Af(ii-1,:);
    a1 = exp(D.*dz/2.0) .* (dz * operator_nieliniowy1(At(ii - 1,:), hR_f, fr, w, t, gamma));
    
    Nl2 = iFT_scal(w,AI + 1/2.0 * a1);
    a2 = dz * operator_nieliniowy1(Nl2, hR_f, fr, w, t, gamma);
   
    Nl3 = iFT_scal(w,AI + 1/2.0 * a2);
    a3 = dz * operator_nieliniowy1(Nl3, hR_f, fr, w, t, gamma);
   
    Nl4 = iFT_scal(w,exp(D.*dz/2.0) .* (AI + a3));
    a4 = dz * operator_nieliniowy1(Nl4, hR_f, fr, w, t, gamma);
   
    Af(ii,:) = exp(D.*dz/2.0) .* (AI + 1/6.0 *(a1 + 2 * a2 + 2 * a3)) + 1/6.0 .* a4;
    At(ii,:) = iFT_scal(w,Af(ii,:));
   
	propagedlength = propagedlength + dz;
    len(ii) = propagedlength;
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /z );
   
	ii = ii + 1;
   
end

tx = toc;

fprintf(1, '\n\nCzas trwania symulacji (s) = ');
fprintf(1, '%5.2f%', tx );

fprintf(1, '\n\nIlosc iteracji = ');
fprintf(1, '%5.0f%', ii );
fprintf(1, '\n\n');

%--------------------------------------------------------------------------
%Wykresy
%--------------------------------------------------------------------------

figure(1)
lIT = 10*log10(abs(At).^2);              % skala logarytmiczna mocy
mlIT = max(max(lIT));                    % maksimum mocy
nlIT= lIT - mlIT;
x = t;
x(t<-0.3) = 0;
x(t>0.3) = 0;
mesh(x, len./(pi/2 * Ld), abs(At).^2, 'FaceLighting','gouraud','LineWidth',0.3);
xlabel('T [ps]'); ylabel('z/z_{sol}'); zlabel('|A(z,T)|^2');
set(gca,'xlim',[-0.3 0.3])
%print('-f1','Dudley_3D_vol2','-dpng')

figure(2)
pcolor(t, len./(pi/2 * Ld), nlIT);           
caxis([-40, 0]);  xlim([-0.3,0.3]); shading interp;
xlabel('T [ps]'); ylabel('z/z_{sol}');
colorbar
%print('-f2','Dudley_time_vol2','-dpng')

figure(3)
w = w + w_c;
lambda = 2*pi*c./w *10^12;                                   %siatka punktów d³ugoœci fali [nm]
lIW = 10*log10(abs(Af).^2);                                  %skala logarytmiczna intensywoœci spektralnej
mlIW = max(max(lIW));                                        %mmaksymalna wartoœæ
nlIW= lIW - mlIW;
pcolor(lambda, len./(pi/2 * Ld), nlIW); 
caxis([-40, 0]); xlim([600 1300])
shading interp; 
xlabel('\lambda [nm]'); ylabel('z/z_{sol}');
colorbar
%print('-f3','Dudley_lambda_vol2','-dpng')


