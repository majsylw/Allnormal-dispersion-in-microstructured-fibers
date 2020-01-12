function f_w = FT_scal(t,f)
%Funkcja wyliczaj�ca transformat� fouriera wraz z przeskalowaniem przekszta�canych wielko�ci
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
% 
% Dane wej�ciowe:
% 
% At    A(z,t) - wektor obwiednii amplitudy w danym punkcie z,
%                dla ca�ej siatki czasowej 
% t     wektor wsp�cz�dnych czasowych
%
% Dane wyj�ciowe:
% 
% Nl    A(z,omega) - wektor tego samego rozmiaru co t, 
%       reprezentuj�cy warto�� obwiedni amplitudy dla danego z w domenie cz�stotliwo�ci

Nt = length (t);            %Liczba punkt�w siatki czasu
dt = abs(t(2) - t(1));      %Krok czasowy

f_w = fftshift(ifft(fftshift(f))) * dt * Nt;
%f= fftshift (f);            %t=0-element
%f_w =ifft(f)*Nt*Dt;         %Transformata i skalowanie
%f_w = ifftshift( f_w );     %omega = 0-element

end