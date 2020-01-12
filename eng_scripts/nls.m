function [v,cz] = nls(d,g,w,af,t,h)
%
%Funkcja wyliczaj�ca pojedy�czy krok algorytmu SSFM do rozwi�zywania NLSE
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
% 
% Dane wej�ciowe:
% 
% af    A(z,omega) - wektor obwiednii amplitudy w danym punkcie z,
%                    dla ca�ej siatki cz�stotliwo�ciowej 
% d  	operator dyspersji
% w    	wektor wsp�cz�dnych cz�stotliwo�ciowych
% t     wektor wsp�cz�dnych czasowych
% h	krok po�o�eniowy
% g	wsp�czynnik nieliniowo�ci

% Dane wyj�ciowe:
% 
% v     A(z+h,omega) - wektor obwiednii amplitudy w punkcie z + h,
%                      dla ca�ej siatki cz�stotliwo�ciowej 
% cz    A(z+h,t) - wektor obwiednii amplitudy w punkcie z + h,
%                  dla ca�ej siatki czasowej 

    v = af .* exp(d .* 0.5 * h);
    v = iFT_scal(w,v);
    v = v .* exp(1i * g * abs(v).^2 * h);
    v = FT_scal(t,v);
    v = v .* exp(d .* 0.5 * h);
    cz = iFT_scal(w,v);

end