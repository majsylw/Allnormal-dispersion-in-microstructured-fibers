function f = iFT_scal(omega,f_w)
%Funkcja wyliczaj¹ca odwrotn¹ transformatê fouriera wraz z przeskalowaniem przekszta³canych wielkoœci
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
% 
% Dane wejœciowe:
% 
% f_w    	A(z,omega) - wektor obwiednii amplitudy w danym punkcie z,
%                        dla ca³ej siatki czêstotliwoœciowej
% omega     wektor wspó³czêdnych czêstotliwoœciowych
%
% Dane wyjœciowe:
% 
% f     A(z,t) - wektor tego samego rozmiaru co omega, 
%       reprezentuj¹cy wartoœæ obwiedni amplitudy dla danego z w domenie czasu

Domega = abs( omega (2) -omega (1));     %Krok czêstotliwoœciowy

f = fftshift(fft(fftshift(f_w))) * Domega/(2* pi);
%f_w = fftshift ( f_w );
%f_w = fft( f_w )* Domega/(2* pi);
%f= ifftshift ( f_w );

end