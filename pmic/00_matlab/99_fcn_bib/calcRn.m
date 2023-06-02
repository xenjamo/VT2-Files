function [Rn] = calcRn(N, alpha)

% Drehung um Richtungsvektor N
% M.E.Peter, 03.10.2012 
%
% N     := [n1_vec n2_vec ... nn_vec] in R^(3 x n) (Richtungsvektor)
% alpha := [alpha1 alpha2 ... alphan]^T (Winkel in rad) in R^(n)
%
% Ist N ein 3x1 Zeilenvektor und alpha ein Skalar, wird Drehmatrix um N
% um den winkel alpha bestummen. Drehrichtung mathematisch positiv.
%
% Die Richtungsvektoren in N werden auf die Länge 1 normiert.
%
% Bsp: Rx = calc_Rn([1 0 0].', alpha) ist Rotation um x-Achse
%
% Ist N eine Matrix und alpha ein Vektor mit entsprechender Dimension wird 
% Rn wie nachfolgend berechnet:
% Rn(n1_vec, alpha(1))*Rn(n2_vec, alpha(2))*...*Rn(n3_vec, alpha(3))
%
% Bsp: Rxyz = calc_Rn(eye,alpha) ist der Rheie nach Rotation um die z-,
% y- und x-Achse (Eulerwinkel)
%
% Referenzes: http://de.wikipedia.org/wiki/Rotationsmatrix

%%

[m_N,n_N] = size(N);
[m_alpha,n_alpha] = size(alpha);

if m_N ~= 3 || n_N ~= m_alpha
    error('N muss Matrix der Dimension R^(3 x dim(alpha)) sein');
end
if n_alpha > 1
    error('alpha muss Zeilenvektor der Dimension n, wobei N in R^(3 x n)');
end

Rn = eye(3);
for k = 1:n_N

    n = N(:,k)/norm(N(:,k));
    c = cos(alpha(k));
    s = sin(alpha(k));

    Rn_dum = [[n(1)^2*(1 - c) + c,         n(1)*n(2)*(1 - c) - n(3)*s, n(1)*n(3)*(1 - c) + n(2)*s];...
              [n(1)*n(2)*(1 - c) + n(3)*s, n(2)^2*(1 - c) + c,         n(2)*n(3)*(1 - c) - n(1)*s];...
              [n(1)*n(3)*(1 - c) - n(2)*s, n(2)*n(3)*(1 - c) + n(1)*s, n(3)^2*(1 - c) + c]];

    Rn = Rn_dum*Rn;

end

return


