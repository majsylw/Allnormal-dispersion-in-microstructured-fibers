%Skrypt wyliczaj¹cy wartoœci sta³ej propagacji w zale¿noœci od d³ugoœci
%fali oraz jej dwóch pierwszych pochodnych na podstawie pliku neff.
%Dodatkowo wyznaczono Aeff(lambda).
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.

clear;
clc;
format long e

load neff;

c = 2.99792458*1e-7;                             % prêdkoœæ œwiat³a w pró¿ni [km/ps]

lambda = neff(:,1)' * 10^(-3);                   % lambda [km]
omega = 2.0 * pi * c./lambda;                    % czêstoœæ [2*pi*THz]

nef = neff(:,2)';                                % efektywny wspó³czynnik za³amania
Aeff = neff(:,3)'*10^(-18);                      % efektywne pole modu [km^2]
Beta = nef.*omega/c; 


for ii = 1:length(Beta)-1
    
        beta1(ii) = (Beta(ii+1)-Beta(ii))/(omega(ii+1)-omega(ii)); 
end

for ii = 1:length(beta1)-1
    
        beta2(ii) = (beta1(ii+1)-beta1(ii))/(omega(ii+1)-omega(ii)); 
end

D = -2*pi*c./lambda(1:length(beta1)-1).^2 .* beta2*10^(-12);



figure(20)
plot(omega,nef,'b');
xlabel('\omega [THz]','FontSize',14);ylabel('n_{eff}','FontSize',14);
%print('-f20','neff1_nefw','-dpng')
figure(21)
plot(lambda*10^12,nef, 'm');
xlabel('\lambda [nm]');ylabel('n_{eff}');
%print('-f21','neff1_nefl','-dpng')
figure(22)
plot(lambda*10^12,Beta,'c');
xlabel('\lambda [THz]');ylabel('\beta [1/km]');
%print('-f22','neff1_betal','-dpng')
figure(23)
plot(lambda(1:length(Beta)-1)*10^12,beta1,'k');
xlabel('\lambda [nm]');ylabel('\beta_1 [ps/km]');
%print('-f23','neff1_B1l','-dpng')
figure(24)
plot(lambda(1:length(beta1)-1)*10^12,beta2,'b');
xlabel('\lambda [nm]');ylabel('\beta_2 [ps^2/km]');
%print('-f24','neff1_B2l','-dpng')
figure(25)
plot(lambda(1:length(beta1)-1)*10^12,D,'r',lambda(1:length(beta1)-1)*10^12,0,'k');
ylim([-200, 0]);
xlim([500 3000])
set(gca,'fontsize',14) 
xlabel('\lambda [nm]','FontSize',20);ylabel('D [ps/km/nm]','FontSize',20);
%print('-f25','neff1_D','-dpng')
figure(26)
plot(lambda*10^12,Aeff*10^18,'b');
ylim([0 200]);
xlim([500 3000])
set(gca,'fontsize',14) 
xlabel('\lambda [nm]','FontSize',20);ylabel('A_{eff} [\mum^2]','FontSize',20);
%print('-f26','neff1_Aeff','-dpng')

