%Symulacja ukazuj�ca proces tworzenia solitonu rz�du 5.
%
% 
% 			      d(A(z,t))/dt = ()D+N)A(z,t)
% 			gdzie N jest operatorem nieliniowo�ci, D operatorem dyspersji
%
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.

clear;
clc;

%--------------------------------------------------------------------------
%Parametry numeryczne
%--------------------------------------------------------------------------

T = 12.5;                                   %szeroko�� okna czasowego [ps]
N = 2^13;                                   %ilo�� punkt�w siatki czasowej
M = 700;                                    %ilo�� punkt�w siatki po�o�eniowej
dt = T/N;                                   %krok czasowy [ps]
dw = 1/T;                                   %krok cz�stotliwo�ciowy [THz]

t = zeros(1,N);                             %siatka punkt�w czasowych
w = zeros(1,N);                             %siatka punkt�w cz�stotliwo�ciowych
for i = 1 : N
    t(i) = -T/2 + i*dt;
    w(i) = 2 * pi * (-N/2 + i) * dw;
end

%--------------------------------------------------------------------------
% Warunki pocz�tkowe - impuls pocz�tkowy
%--------------------------------------------------------------------------

Af = zeros(M,N);                            %A(z,w)
At = zeros(M,N);                            %A(z,t)

P0 = 1250;                                  %maksimum mocy impulsu pocz�tkowego [W]
T0 = 0.0284;                                %czas trwania impulsu [ps]

%Impuls pocz�tkowy w domenie czasowej (sekans hiperboliczny)
At(1,:) = sqrt(P0) .* sech(t./T0);

%Impuls pocz�tkowy w domenie cz�stotliwo�ciowej
Af(1,:) = fftshift(ifft(fftshift(At(1,:)))) * dt * N;

%--------------------------------------------------------------------------
%Parametry �wiat�owodu
%--------------------------------------------------------------------------
Beta = [0 -11.83];                           %pochodne beta(i) [ps^i/km]
m = length(Beta);                            %liczba pochodnych beta
straty = 0;                                  %straty [dB/km]
alpha = log(10.^(straty/10));                %wsp�czynnik t�umienno�ci
gamma = 294;                                 %wsp�czynnik nieliniowo�ci [1/W/km]

Ld = (T0^2)/(abs(Beta(2)));                  %droga dyspersji [km]
Ln = 1/(P0 * gamma);                         %droga nieliniowo�ci [km]

c = 2.99792458 * 1e-7;                       %pr�dko�� �wiat�a w pr�ni [km/ps]
lambda_centralna = 835 * 1e-12;              %centralna d�ugo�� fali [km]
w_c = (2.0*pi*c)/lambda_centralna;           %reference frequency [2*pi*THz]

z = (pi/2) * Ld;                         %d�ugo�� �wiat�owodu [km]
dz = z/(M-1);                                %krok po�o�eniowy [km]
len = 0:dz:z;                                %siatka punkt�w po�o�eniowych

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
x = t./T0;
x(t./T0<-2) = 0;
x(t./T0>2) = 0;
mesh(x, len./z, abs(At).^2./P0, 'FaceLighting','gouraud','LineWidth',0.3);
caxis([0 5])
xlabel('T/T_0','fontsize',16); ylabel('z/z_{sol}','fontsize',16); zlabel('|A(z,t)|^2/P_0','fontsize',16);
set(gca,'xlim',[-2 2])
set(gca,'fontsize',14);
print('-f1','soliton5th_3D2','-dpng')