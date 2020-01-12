function Nl = operator_nieliniowy2(Af, hR_f, fr, v, t, v0,gamma) 
 
%Funkcja wyliczaj¹ca wartoœæ operatora nieliniowoœci N(A(z,omega)) = N(IFFT(Af)) 
%w domenie czêstotliwoœci
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
%
%   		Nl(At) = i * gamma * (1 + v/v0)* FFT(At((1-fr)*abs(At)^2 +
%                  + fr * IFFT(hr_f * FFT(abs(At)^2))))
% 
% Dane wejœciowe:
% 
% Af    A(z,omega) - wektor obwiedni amplitudy w danym punkcie z,
%                    dla ca³ej siatki czêstotliwoœciowej 
% hR_f  wektor wartoœci funkcji 
%       hR(t>=0) = (tau1^2 + tau2^2)/(tau1*tau2^2)*exp(-t/tau2)*sin(t/tau1)
%		w domenie czêstotliwoœci
% fr  	maksimum wzmocnienia Ramanowskiego
% v    	wektor wspó³czêdnych czasowych
% t     wektor wspó³czêdnych czasowych
% v0	centralna czêstotliwoœæ impulsu
% gamma	wspó³czynnik nieliniowoœci

% Dane wyjœciowe:
% 
% Nl    wektor tego samego rozmiaru co t, 
%       reprezentuj¹cy wartoœæ N(A(z,omega)) = N(IFFT(Af))

    At = iFT_scal(v,Af);
    A2t = abs(At).^2;
    A2w = FT_scal(t,A2t);
    
    A2Rw = hR_f .* A2w;
    A2Rt = iFT_scal(v,A2Rw);
    
    A3t = At .* ((1 - fr) * A2t + fr * A2Rt);
    A3w = FT_scal(t,A3t);
    
    Nl = 1i * gamma .* (1 + v/v0) .* A3w;
    %Nl = 1i * gamma .* (1 + 0) .* A3w; %w przypadku nieuwzglêdniania samostromoœci
    
    nan_region = find(isnan(Nl));
    Nl(nan_region) = 1e-100;
	
end
