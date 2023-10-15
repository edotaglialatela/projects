function [F] = forcefunct(answer, F0, l_2, t, omega)
if answer==1 %constant
    F = F0.*[1;l_2;0;0];
elseif answer==2 %alternating
    F = F0.*[1;l_2;0;0].*sin(omega.*t);
elseif answer==3 %start
    F=F0.*[1;l_2;0;0].*sin(omega.*t).*(1-exp(-(100*t^2)));
end
end