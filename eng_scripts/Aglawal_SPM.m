
%Symulacja ukazuj¹ca proces samomodulacji fazy
% Agrawal, G. P., Nonlinear Fiber Optics
% 
% 			      d(A(z,t))/dz = NA(z,t)
% 			gdzie N jest operatorem nieliniowoœci
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.


clear;
clc;

%--------------------------------------------------------------------------
%Parametry numeryczne
%--------------------------------------------------------------------------

T = 48000;                  %szerokoœæ okna czasowego [ps]
N = 2^13;                   %iloœæ punktów siatki czasowej
M = 50;                     %iloœæ punktów siatki po³o¿eniowej
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

Af = zeros(M,N);            %A(z,w)
At = zeros(M,N);            %A(z,t)

P0 = 1;                     %maksimum mocy impulsu pocz¹tkowego [W]
T0 = 200;                   %czas trwania impulsu [ps]

%Impuls poczatkowy w domenie czasowej (gauss)
At(1,:) = sqrt(P0) .* exp(-1/2 .* (t./T0) .^ 2);
%Impuls poczatkowy w domenie czêstotliwoœciowej
Af(1,:) = fftshift(ifft(fftshift(At(1,:)))) * dt * N;
%Af(1,:) = FT_scal(t,At(1,:));

%--------------------------------------------------------------------------
%Parametry œwiat³owodu
%--------------------------------------------------------------------------

gamma = 3;                  %wspó³czynnik nieliniowoœci [1/W/km]
Ln = 1/(P0 * gamma);        %droga nieliniowoœci [km]
z = 3.5 * pi * Ln;          %d³ugoœæ œwiat³owodu [km]
dz = z/(M-1);               %krok po³o¿eniowy [km]
len = 0:dz:z;               %siatka punktów po³o¿eniowych

%--------------------------------------------------------------------------

for ii = 2:M
    At(ii, :) = At(ii-1,:) .* exp(1i * gamma * abs(At(ii-1,:)).^2 * dz);
    %Af(ii,:) = FT_scal(t,At(ii,:));
    Af(ii,:) = fftshift(ifft(fftshift(At(ii,:)))) * dt * N;
end


%--------------------------------------------------------------------------
%Wykresy
%--------------------------------------------------------------------------

figure(1)
x = w./(2 * pi) * T0;
x(w./(2 * pi) * T0<-3) = 0;
x(w./(2 * pi) * T0>3) = 0;
mesh(x, len./(Ln * pi), abs(Af).^2./max(abs(Af(1,:)).^2), 'FaceLighting','gouraud','LineWidth',0.3);
xlabel('\omega/(2\pi) * T_0'); ylabel('z/(\piL_N)'); zlabel('Spektralna gêstoœæ energii');
axis([-3 3 0 3.5 0 1]);

figure(2);
contourf( w./(2 * pi) * T0, len./(Ln * pi), abs(Af).^2./max(abs(Af(1,:)).^2), 100, 'LineStyle', 'none')
caxis([0 1])
xlim([-3 3])
xlabel('\omega/(2\pi) * T_0'); ylabel('z/(\piL_N)'); zlabel('Spektralna gêstoœæ energii');
colorbar

figure(3);
plot(w./(2 * pi) * T0, abs(Af(1, : )).^2/max(abs(Af(1,:)).^2), 'red', w./(2 * pi) * T0, abs(Af(10,:)).^2/max(abs(Af(1,:)).^2), 'blue', w./(2 * pi) * T0, abs(Af(20,:)).^2/max(abs(Af(1,:)).^2), 'green', w./(2 * pi) * T0, abs(Af(30,:)).^2/max(abs(Af(1,:)).^2), 'k', w./(2 * pi) * T0, abs(Af(40,:)).^2/max(abs(Af(1,:)).^2), 'yellow', w./(2 * pi) * T0, abs(Af(50,:)).^2/max(abs(Af(1,:)).^2), 'c');
xlim([-2 2])
xlabel('\omega/(2\pi) * T_0','fontsize',16); ylabel('Spektralna gêstoœæ energii','fontsize',16);
set(gca,'fontsize',14);
%print('-f3','3','-dpng')

