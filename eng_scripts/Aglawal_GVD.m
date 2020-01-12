
%Symulacja ukazuj¹ca wp³ywu dyspersji na poszerzenie impulsu
% Agrawal, G. P., Nonlinear Fiber Optics.
% 
% 		   d(A(z,omega))/dz = D(omega)A(z,omega)
% 			gdzie D jest operatorem dyspersji
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.

clear;
clc;
%--------------------------------------------------------------------------
%Parametry numeryczne
%--------------------------------------------------------------------------

T = 48000;                      %szerokoœæ okna czasowego [ps]
N = 2^13;                       %iloœæ punktów siatki czasowej
M = 40;                         %iloœæ punktów siatki po³o¿eniowej
dt = T/N;                       %krok czasowy [ps]
dw = 1/T;                       %krok czêstotliwoœciowy [THz]
    
t = zeros(1,N);                 %siatka punktów czasowych
w = zeros(1,N);                 %siatka punktów czêstotliwoœciowych
for i = 1 : N
    t(i) = -T/2 + i*dt;
    w(i) = 2 * pi * (-N/2 + i) * dw;
end

%--------------------------------------------------------------------------
% Warunki pocz¹tkowe - impuls pocz¹tkowy
%--------------------------------------------------------------------------

Af = zeros(M,N); %A(z,w)
At = zeros(M,N); %A(z,t)

P0 = 1;                         %maksimum mocy impulsu pocz¹tkowego [W]
T0 = 250;                       %czas trwania impulsu [ps]

%Impuls poczatkowy w domenie czasowej (gauss)
At(1,:) = sqrt(P0) .* exp(-1/2 .* (t./T0) .^ 2);

%Impuls poczatkowy w domenie czêstotliwoœciowej
Af(1,:) = fftshift(ifft(fftshift(At(1,:)))) * dt * N;
%Af(1,:) = FT_scal(t,At(1,:));

%--------------------------------------------------------------------------
%Parametry œwiat³owodu
%--------------------------------------------------------------------------
Beta = [0 20 0 0];              %pochodne beta(i) [ps^i/km]
m = length(Beta);               %liczba pochodnych beta
straty = 0;                     %straty [dB/km]
alpha = log(10.^(straty/10));   %wspó³czynnik t³umiennoœci            
Ld = (T0^2)/(abs(Beta(2)));     %droga dyspersji [km]
z = 4 * Ld;                     %d³ugoœæ œwiat³owodu [km]
dz = z/(M-1);                   %krok po³o¿eniowy [km]
len = 0:dz:z;                   %siatka punktów po³o¿eniowych

%--------------------------------------------------------------------------
%Operator dyspersji
%--------------------------------------------------------------------------

B = 0;                          %Beta(w)

for k = 2:m
    B = B + Beta(k)/factorial(k) * w.^k;
end

D = 1i * B - alpha/2;           %operator liniowy

%--------------------------------------------------------------------------

for ii = 2:M
    Af(ii,:) = Af(ii-1,:) .* exp(D .* dz);
    At(ii,:) = fftshift(fft(fftshift(Af(ii,:)))) * dw;
    %At( ii, : ) = iFT_scal(w,Af( ii, : ));
end

%--------------------------------------------------------------------------
%Wykresy
%--------------------------------------------------------------------------

figure( 1 )
x = t./T0;
x(t./T0<-10) = 0;
x(t./T0>10) = 0;
mesh(x, len./Ld, abs(At).^2./P0, 'FaceLighting','gouraud','LineWidth',0.3);
xlabel('t/T_0'); ylabel('z/L_D'); zlabel('|A(z,T)|^2');
set(gca,'xlim',[-10, 10],'ylim',[0, 4],'zlim',[0, 1])

figure( 2 )
contourf( t./T0, len./Ld, abs( At ).^2./P0, 100, 'LineStyle', 'none' )
caxis([0 1])
xlim([-10 10])
xlabel('T/T_0'); ylabel('z/L_D');
colorbar

figure( 3 )
plot(t./T0, abs(At(1,:)).^2/P0, 'red', t./T0, abs(At(20,:)).^2/P0, 'blue', t./T0, abs(At(40,:)).^2/P0, 'green')
xlim([-10 10])
xlabel('T/T_0'); ylabel('|A(z,T)|^2');
