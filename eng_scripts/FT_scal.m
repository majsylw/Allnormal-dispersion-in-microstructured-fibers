function f_w = FT_scal(t,f)
%Funkcja wyliczaj¹ca transformatê fouriera wraz z przeskalowaniem przekszta³canych wielkoœci
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
% 
% Dane wejœciowe:
% 
% At    A(z,t) - wektor obwiednii amplitudy w danym punkcie z,
%                dla ca³ej siatki czasowej 
% t     wektor wspó³czêdnych czasowych
%
% Dane wyjœciowe:
% 
% Nl    A(z,omega) - wektor tego samego rozmiaru co t, 
%       reprezentuj¹cy wartoœæ obwiedni amplitudy dla danego z w domenie czêstotliwoœci

Nt = length (t);            %Liczba punktów siatki czasu
dt = abs(t(2) - t(1));      %Krok czasowy

f_w = fftshift(ifft(fftshift(f))) * dt * Nt;
%f= fftshift (f);            %t=0-element
%f_w =ifft(f)*Nt*Dt;         %Transformata i skalowanie
%f_w = ifftshift( f_w );     %omega = 0-element

end