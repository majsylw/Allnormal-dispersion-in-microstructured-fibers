function f = iFT_scal(omega,f_w)
%Funkcja wyliczaj�ca odwrotn� transformat� fouriera wraz z przeskalowaniem przekszta�canych wielko�ci
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
% 
% Dane wej�ciowe:
% 
% f_w    	A(z,omega) - wektor obwiednii amplitudy w danym punkcie z,
%                        dla ca�ej siatki cz�stotliwo�ciowej
% omega     wektor wsp�cz�dnych cz�stotliwo�ciowych
%
% Dane wyj�ciowe:
% 
% f     A(z,t) - wektor tego samego rozmiaru co omega, 
%       reprezentuj�cy warto�� obwiedni amplitudy dla danego z w domenie czasu

Domega = abs( omega (2) -omega (1));     %Krok cz�stotliwo�ciowy

f = fftshift(fft(fftshift(f_w))) * Domega/(2* pi);
%f_w = fftshift ( f_w );
%f_w = fft( f_w )* Domega/(2* pi);
%f= ifftshift ( f_w );

end