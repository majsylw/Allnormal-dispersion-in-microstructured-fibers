%Symulacja ukazuj¹ca proces tworzenia solitonu rzêdu 3.
%Dudley, J. M., G. Genty oraz S. Coen, 
%„Supercontinuum generation in photoniccrystal fiber,” 
%Review of Modern Physcis, 78(4), str. 1135–1184, pazdziernik 2006.
% 
% NLSE
% 	   d(A(z,omega))/dt = (D(omega)+N(A(z,t)))A(z,omega)
% 	gdzie N jest operatorem nieliniowoœci, D operatorem dyspersji
%
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.

clear;
clc;

%--------------------------------------------------------------------------
%Parametry numeryczne
%--------------------------------------------------------------------------

T = 12.5;                                   %szerokoœæ okna czasowego [ps]
N = 2^13;                                   %iloœæ punktów siatki czasowej
M = 500;                                    %iloœæ punktów siatki po³o¿eniowej
dt = T/N;                                   %krok czasowy [ps]
dw = 1/T;                                   %krok czêstotliwoœciowy [THz]

t = zeros(1,N);                             %siatka punktów czasowych
w = zeros(1,N);                             %siatka punktów czêstotliwoœciowych
for i = 1 : N
    t(i) = -T/2 + i*dt;
    w(i) = 2 * pi * (-N/2 + i) * dw;
end

%--------------------------------------------------------------------------
% Warunki pocz¹tkowe - impuls pocz¹tkowy
%--------------------------------------------------------------------------

Af = zeros(M,N);                            %A(z,w)
At = zeros(M,N);                            %A(z,t)

P0 = 1250;                                  %maksimum mocy impulsu pocz¹tkowego [W]
T0 = 0.0284;                                %czas trwania impulsu [ps]

%Impuls pocz¹tkowy w domenie czasowej (sekans hiperboliczny)
At(1,:) = sqrt(P0) .* sech(t./T0);

%Impuls pocz¹tkowy w domenie czêstotliwoœciowej
Af(1,:) = fftshift(ifft(fftshift(At(1,:)))) * dt * N;

%--------------------------------------------------------------------------
%Parametry œwiat³owodu
%--------------------------------------------------------------------------
Beta = [0 -11.83];                           %pochodne beta(i) [ps^i/km]
m = length(Beta);                            %liczba pochodnych beta
straty = 0;                                  %straty [dB/km]
alpha = log(10.^(straty/10));                %wspó³czynnik t³umiennoœci
gamma = 110;                                 %wspó³czynnik nieliniowoœci [1/W/km]

Ld = (T0^2)/(abs(Beta(2)));                  %droga dyspersji [km]
Ln = 1/(P0 * gamma);                         %droga nieliniowoœci [km]

c = 2.99792458 * 1e-7;                       %prêdkoœæ œwiat³a w pró¿ni [km/ps]
lambda_centralna = 835 * 1e-12;              %centralna d³ugoœæ fali [km]
w_c = (2.0*pi*c)/lambda_centralna;           %czêstoœæ centralna [2*pi*THz]

z = 2 * (pi/2) * Ld;                         %d³ugoœæ œwiat³owodu [km]
dz = z/(M-1);                                %krok po³o¿eniowy [km]
len = 0:dz:z;                                %siatka punktów po³o¿eniowych

%--------------------------------------------------------------------------
%Operator dyspersji
%--------------------------------------------------------------------------

B = 0;                                      %Beta(w)

for k = 2:m
    B = B + Beta(k)/factorial(k) * w.^k;
end

D = 1i * B - alpha/2;                       %operator liniowy


%--------------------------------------------------------------------------

for ii = 2:M
    [Af(ii,:),At(ii,:)] = nls(D,gamma,w,Af(ii - 1,:),t,dz);    
end

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

figure(2)
pcolor(t, len./(pi/2 * Ld), nlIT);           
caxis([-40, 0]);  xlim([-0.3,0.3]); shading interp;
xlabel('T [ps]'); ylabel('z/z_{sol}');
colorbar

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
