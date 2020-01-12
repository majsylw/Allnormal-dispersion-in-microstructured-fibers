function Nl = operator_nieliniowy2(Af, hR_f, fr, v, t, v0,gamma) 
 
%Funkcja wyliczająca wartość operatora nieliniowości N(A(z,omega)) = N(IFFT(Af)) 
%w domenie częstotliwości
%
% Autor: Sylwia Majchrowska
% 10 grudnia 2015r.
%
%   		Nl(At) = i * gamma * (1 + v/v0)* FFT(At((1-fr)*abs(At)^2 +
%                  + fr * IFFT(hr_f * FFT(abs(At)^2))))
% 
% Dane wejściowe:
% 
% Af    A(z,omega) - wektor obwiedni amplitudy w danym punkcie z,
%                    dla całej siatki częstotliwościowej 
% hR_f  wektor wartości funkcji 
%       hR(t>=0) = (tau1^2 + tau2^2)/(tau1*tau2^2)*exp(-t/tau2)*sin(t/tau1)
%		w domenie częstotliwości
% fr  	maksimum wzmocnienia Ramanowskiego
% v    	wektor współczędnych czasowych
% t     wektor współczędnych czasowych
% v0	centralna częstotliwość impulsu
% gamma	współczynnik nieliniowości

% Dane wyjściowe:
% 
% Nl    wektor tego samego rozmiaru co t, 
%       reprezentujący wartość N(A(z,omega)) = N(IFFT(Af))

    At = iFT_scal(v,Af);
    A2t = abs(At).^2;
    A2w = FT_scal(t,A2t);
    
    A2Rw = hR_f .* A2w;
    A2Rt = iFT_scal(v,A2Rw);
    
    A3t = At .* ((1 - fr) * A2t + fr * A2Rt);
    A3w = FT_scal(t,A3t);
    
    Nl = 1i * gamma .* (1 + v/v0) .* A3w;
    %Nl = 1i * gamma .* (1 + 0) .* A3w; %w przypadku nieuwzględniania samostromości
    
    nan_region = find(isnan(Nl));
    Nl(nan_region) = 1e-100;
	
end
