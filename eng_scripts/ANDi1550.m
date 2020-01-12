%Generacja supercontinuum w re¿imie dyspersji ca³kowitonormalnej
%K. L. Tarnowski, W. Urbanczyk,
%„All-normal dispersion hole-assisted silica fibers for
%generation of supercontinuum reaching mid-infrared”, 
%Optics Express, w przygotowaniu.
%
% 
% 	  d(A(z,omega))/dz = (N(A(z,omega))+D(omega))A(z,omega)
% 			gdzie N jest operatorem nieliniowoœci, a D dyspersji
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.

clear;
clc;
format long e

c = 2.99792458*1e-7;                         %prêdkoœæ œwiat³a w pró¿ni [km/ps]
lambda_centralna = 1550*1e-12;               %centralna d³ugoœæ fali [km]
w_c = (2.0*pi*c)/lambda_centralna;           %centralna czêstoœæ [2*pi*THz]

%--------------------------------------------------------------------------
%Parametry numeryczne
%--------------------------------------------------------------------------

T = 22.6;                   %szerokoœæ okna czasowego [ps]
N = 2^13;                   %liczba punktów siatki czasowej               
dt = T/N;                   %krok czasowy [ps]
dw = 1/T;                   %krok czêstotliwoœciowy [THz]

t = linspace(-T/2, T/2, N); %siatka punktów czasowych
w = 2*pi*(-N/2:N/2-1) * dw; %siatka punktów czêstotliwoœciowych

%--------------------------------------------------------------------------
% Warunki pocz¹tkowe - impuls pocz¹tkowy
%--------------------------------------------------------------------------

%Af = A(z,w)
%At = A(z,t)

P0 = 10000;                                     % maksimum mocy impulsu pocz¹tkowego [W]
TFWHM = 0.065;                                  % szerokoœæ po³ówkowa [ps]
T0 = TFWHM/(2*log(1+sqrt(2)));                  % czas trwania impulsu [ps]

%Impuls pocz¹tkowy w domenie czasowej (sekans hiperboliczny)
At(1,:) = sqrt(P0) .* sech(t./T0);

%Impuls pocz¹tkowy w domenie czêstotliwoœciowej
Af(1,:) = FT_scal(t,At(1,:));

%--------------------------------------------------------------------------
%Parametry œwiat³owodu
%--------------------------------------------------------------------------
alpha = 0;                                       % wspó³czynnik t³umiennoœci [W/km]

fr = 0.18;                                       
tau1 = 0.0122;                                   % [ps]
tau2 = 0.032;                                    % [ps]

z = 10^(-3);                                     % d³ugoœæ œwiat³owodu [km]
dz = 2e-7;                                       % krok przestrzenny
%len                                   			 % siatka punktów przestrzennych

%--------------------------------------------------------------------------
%Operator dyspersji
%--------------------------------------------------------------------------
load neff;

lambda = neff(:,1)' * 10^(-3);                   % lambda [km]
omega = 2.0 * pi * c./lambda;                    % czêstoœæ [2*pi*THz]

nef = neff(:,2)';                                % efektywny wspó³czynnik za³amania
Aeff = neff(:,3)'*10^(-18);                      % efektywne pole modu [km^2]
Betta = nef.*omega/c; 

for ii = 1:length(Betta)-1
    
        beta1(ii) = (Betta(ii+1)-Betta(ii))/(omega(ii+1)-omega(ii)); 
end

Beta_0 = Betta(lambda == lambda_centralna);      % B(omega_0)
nef1_0 = nef(lambda == lambda_centralna);        % n_eff(omega_0)
Beta1_0 = beta1(lambda == lambda_centralna);     % B_1(omega_0)
Aeff_0 = Aeff(lambda == lambda_centralna);       % A_eff(omega_0)

om = w + w_c;
llam = 2*pi*c./om * 10^12;                       % siatka punktów d³ugoœci fali [nm]

nef1 = interp1(omega, nef, om);
nan_region = find(isnan(nef1));
nef1(nan_region) = nef1(length(nan_region)+1);

Beta = interp1(omega, Betta, om);
nan_region = find(isnan(Beta));
%Beta(nan_region) = Beta(length(nan_region)+1);
Beta(nan_region) = 40/(1i);

D = - alpha/2 + 1i * (Beta - Beta_0 - Beta1_0 * w);       %operator liniowy

%--------------------------------------------------------------------------
%Wspó³czynnik nieliniowoœci
%--------------------------------------------------------------------------

Aef1 = interp1(omega, Aeff, om);
nan_region = find(isnan(Aef1));
Aef1(nan_region) = 1;

poprawka = (Aef1/Aeff_0).^(-1/4);
gamma = 2.7e-26 * nef1_0 * w_c./(c * nef1 .* sqrt(Aeff_0 * Aef1));    %wspó³czynnik nieliniowoœci [1/W/km]

%poprawka = 1;
%gamma = 2.7e-26 * w_c / (c * Aeff_0);

%--------------------------------------------------------------------------
%Rozpraszanie Ramana
%--------------------------------------------------------------------------

hR_t = (tau1^2+tau2^2)/(tau1*tau2^2)*(exp(-t/tau2).*sin(t/tau1));
hR_t(t<0) = 0;
hR_f = FT_scal(t,hR_t);

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
    
    cc = poprawka .* Af(ii-1,:);
    
    AI = exp(D.*dz/2) .* cc;
    a1 = exp(D.*dz/2) * dz .* operator_nieliniowy2(cc, hR_f, fr, w, t, w_c, gamma);
   
    Nl2 = (AI + 1/2 * a1);
    a2 = dz * operator_nieliniowy2(Nl2, hR_f, fr, w, t, w_c, gamma);
   
    Nl3 = (AI + 1/2 * a2);
    a3 = dz * operator_nieliniowy2(Nl3, hR_f, fr, w, t, w_c, gamma);
   
    Nl4 = (exp(D.*dz/2) .* (AI + a3));
    a4 = dz * operator_nieliniowy2(Nl4, hR_f, fr, w, t, w_c, gamma);
   
    Af(ii,:) = exp(D.*dz/2) .* (AI + 1/6 *(a1 + 2 * a2 + 2 * a3)) + 1/6 .* a4;
    Af(ii,:) = Af(ii,:)./poprawka;
    At(ii,:) = iFT_scal(w,Af(ii,:));
   
	propagedlength = propagedlength + dz;
    len(ii) = propagedlength;
    fprintf(1, '\b\b\b\b\b\b%5.2f%%', propagedlength * 100.0 /z );
   
	ii = ii + 1;
   
end

tx = toc;

fprintf(1, '\n\nCzas trwania symulacji (s) = ');
fprintf(1, '%5.2f%', tx );

fprintf(1, '\n\nLiczba iteracji = ');
fprintf(1, '%5.0f%', ii - 1 );
fprintf(1, '\n\n');
%--------------------------------------------------------------------------
%Wykresy
%--------------------------------------------------------------------------
%save('Af_2e-7_215sup_100kW_65fs.mat', 'Af');
%save('At_1e-6_215sup_10kW_65fs.mat', 'At');
return

figure(1)
lIT = 10*log10(abs(At).^2);              % skala logarytmiczna mocy
mlIT = max(max(lIT));                    % maksimum mocy
nlIT= lIT - mlIT;
x = t;
x(t<-7) = 0;
x(t>7) = 0;
mesh(x, len.*10^3, abs(At).^2, 'FaceLighting','gouraud','LineWidth',0.3);
xlabel('T [ps]','fontsize',16); ylabel('z [m]','fontsize',16); zlabel('|A(z,T)|^2','fontsize',16);
set(gca,'xlim',[-7 7])
set(gcf,'renderer','zbuffer');
set(gca,'fontsize',14);
%print('-f1','ANDI1550_3D_while_FD_1e-6_213sup_10kW_65fse','-dpng')


figure(2)
pcolor(t, len.*10^3, nlIT);           
caxis([-30, 0]); 
xlim([-7,7]); 
shading interp;
xlabel('T [ps]','fontsize',16); ylabel('z [m]','fontsize',16);
colorbar
set(gcf,'renderer','zbuffer');
set(gca,'fontsize',14);
%print('-f2','ANDI1550_time_while_FD_1e-6_213sup_10kW_65fse','-dpng')


figure(3)                                  
lIW = 10*log10(abs(Af).^2);                                  %skala logarytmiczna intensywoœci spektralnej
mlIW = max(max(lIW));                                        %mmaksymalna wartoœæ
nlIW= lIW - mlIW;
pcolor(llam, len.* 10^3, nlIW); 
caxis([-30, 0]); xlim([800 2200])
shading interp; 
xlabel('\lambda [nm]','fontsize',16); ylabel('z [m]','fontsize',16);
colorbar
set(gcf,'renderer','zbuffer');
set(gca,'fontsize',14);
%print('-f3','ANDI1550_lambda_while_FD_1e-6_213sup_10kW_65fse','-dpng')

%clear all;
%close all;
%clc;