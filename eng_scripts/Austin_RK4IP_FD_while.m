%Generacja supercontinuum na podstawie parametrów z 
%Austin, D. R., C. M. de Sterke oraz B. J. Eggleton, 
%„Dispersive wave blueshift in supercontinuum generation”, 
%Optics Express, 14(25), str. 11 997–12 007, Grudzien 2006.
%
% 
% 	  d(A(z,omega))/dz = (N(A(z,omega))+D(omega))A(z,omega)
% 			gdzie N jest operatorem nieliniowoœci, a D operatorem dyspersji
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.

clear;
clc;
format long e

%--------------------------------------------------------------------------
%Parametry numeryczne
%--------------------------------------------------------------------------

T = 12;                     %szerokoœæ okna czasowego [ps]
N = 2^13;                   %iloœæ punktów siatki czasowej
dt = T/N;                   %krok czasowy [ps]
dw = 1/T;                   %krok czêstotliwoœciowy [THz]

t = linspace(-T/2, T/2, N); %siatka punktów czasowych
w = 2*pi*(-N/2:N/2-1) * dw; %siatka punktów czêstotliwoœciowych

%--------------------------------------------------------------------------
% Warunki pocz¹tkowe - impuls pocz¹tkowy
%--------------------------------------------------------------------------

%Af = A(z,w)
%At = A(z,t)

P0 = 1365;                  %maksimum mocy impulsu pocz¹tkowego [W]
T0 = 0.045;                  %czas trwania impulsu [ps]

%Impuls pocz¹tkowy w domenie czasowej (sekans hiperboliczny)
At(1,:) = sqrt(P0) .* sech(t./T0);

%Impuls pocz¹tkowy w domenie czêstotliwoœciowej
Af(1,:) = FT_scal(t,At(1,:));


%--------------------------------------------------------------------------
%Parametry œwiat³owodu
%--------------------------------------------------------------------------
Beta = [0 -10.64 0.04875 8.164e-5];          %pochodne beta(i) [ps^i/km]
m = length(Beta);                            %liczba pochodnych beta
alpha = 0;                                   %wspó³czynnik t³umiennoœci [W/km]
gamma = 96.2;                                %wspó³czynnik nieliniowoœci [1/W/km]

fr = 0.18;                                   %
tau1 = 0.0122;                               %[ps]
tau2 = 0.032;                                %[ps]

c = 2.99792458*1e-7;                         %prêdkoœæ œwiat³a w pró¿ni [km/ps]
lambda_centralna = 810*1e-12;                %centralna d³ugoœæ fali [km]
w_c = (2.0*pi*c)/lambda_centralna;           %reference frequency [2*pi*THz]

Ld = (T0^2)/(abs(Beta(2)));                  %droga dyspersji [km]
Ln = 1/(P0 * gamma);                         %droga nieliniowoœci [km]

z = 10^-4;                                   %d³ugoœæ œwiat³owodu [km]
dz = z/500;                                  %krok przestrzenny
%len = 0:dz:z;                               %siatka punktów przestrzennych

%--------------------------------------------------------------------------
%Operator dyspersji
%--------------------------------------------------------------------------

D = - alpha/2;       %operator liniowy
for k = 2:m
    D = D + 1i * Beta(k) * w.^k/factorial(k);
end

%--------------------------------------------------------------------------
%Rozpraszanie Ramana
%--------------------------------------------------------------------------

hR_t = (tau1^2+tau2^2)/(tau1*tau2^2)*(exp(-t/tau2).*sin(t/tau1));
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
   
    AI = exp(D.*dz/2) .* Af(ii-1,:);
    a1 = exp(D.*dz/2) * dz .* operator_nieliniowy(At(ii - 1,:), hR_f, fr, w, t, w_c, gamma);
   
    Nl2 = iFT_scal(w,AI + 1/2 * a1);
    a2 = dz * operator_nieliniowy(Nl2, hR_f, fr, w, t, w_c, gamma);
   
    Nl3 = iFT_scal(w,AI + 1/2 * a2);
    a3 = dz * operator_nieliniowy(Nl3, hR_f, fr, w, t, w_c, gamma);
   
    Nl4 = iFT_scal(w,exp(D.*dz/2) .* (AI + a3));
    a4 = dz * operator_nieliniowy(Nl4, hR_f, fr, w, t, w_c, gamma);
   
    Af(ii,:) = exp(D.*dz/2) .* (AI + 1/6 *(a1 + 2 * a2 + 2 * a3)) + 1/6 .* a4;
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
fprintf(1, '%5.0f%', ii - 1 );
fprintf(1, '\n\n');
%--------------------------------------------------------------------------
%Wykresy
%--------------------------------------------------------------------------

figure(1)
lIT = 10*log10(abs(At).^2);              % skala logarytmiczna mocy
mlIT = max(max(lIT));                    % maksimum mocy
nlIT= lIT - mlIT;
x = t;
x(t<-0.2) = 0;
x(t>0.5) = 0;
mesh(x, len.*10^3, abs(At).^2, 'FaceLighting','gouraud','LineWidth',0.3);
xlabel('T [ps]'); ylabel('z [m]'); zlabel('|A(z,T)|^2');
set(gca,'xlim',[-0.2 0.5])
set(gcf,'renderer','zbuffer');
%print('-f1','Austin_3D_while_FD_1001','-dpng')

figure(2)
pcolor(t, len.*10^3, nlIT);           
caxis([-40, 0]); 
xlim([-0.2,0.7]); 
shading interp;
xlabel('T [ps]'); ylabel('z [m]');
colorbar
set(gcf,'renderer','zbuffer');
print('-f2','Austin_time_while_FD_1001','-dpng')

figure(3)
w = w + w_c;
lambda = 2*pi*c./w *10^12;                                   %siatka punktów d³ugoœci fali [nm]
lIW = 10*log10(abs(Af).^2);                                  %skala logarytmiczna intensywoœci spektralnej
mlIW = max(max(lIW));                                        %mmaksymalna wartoœæ
nlIW= lIW - mlIW;
pcolor(lambda, len.* 10^3, nlIW); 
caxis([-40, 0]); xlim([600 1000])
shading interp; 
xlabel('\lambda [nm]'); ylabel('z [m]');
colorbar
set(gcf,'renderer','zbuffer');
print('-f3','Austin_lambda_while_FD_1001','-dpng')

clear all;
close all;
clc;