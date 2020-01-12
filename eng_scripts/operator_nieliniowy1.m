function Nl = operator_nieliniowy1(At, hR_f, fr, v, t, gamma) 
 
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
% At    A(z,t) - wektor obwiedni amplitudy w danym punkcie z,
%                dla ca³ej siatki czasowej 
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
    
    dt = t(2) - t(1);
  
    A2t = abs(At).^2;
    A2w = FT_scal(t,A2t);
    
    A2Rw = hR_f .* A2w;
    A2Rt = iFT_scal(v,A2Rw);
    
    A3t = At .* ((1 - fr) * A2t + fr * A2Rt);
    A3w = FT_scal(t,A3t);
    
    Nl = 1i * gamma .* A3w;
    
    nan_region = find(isnan(Nl));
    Nl(nan_region) = 0;
	
end
