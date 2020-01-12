function [v,cz] = nls(d,g,w,af,t,h)
%
%Funkcja wyliczaj¹ca pojedyñczy krok algorytmu SSFM do rozwi¹zywania NLSE
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
% 
% Dane wejœciowe:
% 
% af    A(z,omega) - wektor obwiednii amplitudy w danym punkcie z,
%                    dla ca³ej siatki czêstotliwoœciowej 
% d  	operator dyspersji
% w    	wektor wspó³czêdnych czêstotliwoœciowych
% t     wektor wspó³czêdnych czasowych
% h	krok po³o¿eniowy
% g	wspó³czynnik nieliniowoœci

% Dane wyjœciowe:
% 
% v     A(z+h,omega) - wektor obwiednii amplitudy w punkcie z + h,
%                      dla ca³ej siatki czêstotliwoœciowej 
% cz    A(z+h,t) - wektor obwiednii amplitudy w punkcie z + h,
%                  dla ca³ej siatki czasowej 

    v = af .* exp(d .* 0.5 * h);
    v = iFT_scal(w,v);
    v = v .* exp(1i * g * abs(v).^2 * h);
    v = FT_scal(t,v);
    v = v .* exp(d .* 0.5 * h);
    cz = iFT_scal(w,v);

end