function [H,H2,paras] = PlaneWaveBC(x,i,t,paras)
% global c_eta_0


mag = paras{1};
phi = paras{2};
omega = paras{3};
omega2 =3.738e14*2*pi;
omega3 = 2.738e14*2*pi;
betap = paras{4};
t0 = paras{5};
st = paras{6};
s = paras{7};
x0 = paras{8};
stx = paras{9};
if length(paras) == 10
    type = paras{10};
else 
    type = 's';
end

if length(paras) == 11
    n = paras{11};
else 
    n = 2;
end


if type == 's'
    if s && t > t0 || st == 0
        H = mag.*sin(omega*t + phi + x*betap);
        H2 = mag.*sin(omega2*t + phi + x*betap);
    elseif st < 0 && t > t0
        H = mag.*sin(omega*t + phi + x*betap);
        H2 = mag.*sin(omega2*t + phi + x*betap);
    else
        H = mag.*sin(omega*t + phi + x*betap)*exp(-((t-t0)/st)^n);
        H2 = mag.*sin(omega2*t + phi + x*betap)*exp(-((t-t0)/st)^n);
    end
elseif type == 'p'
    
    if t > t0 && t < t0+st
        H = mag;
    else 
        H = 0;
    end
    
end

    if stx ~= 0
        H = H.*exp(-((x-x0)/stx).^2);
        H2 = H2.*exp(-((x-x0)/stx).^2);
    end
% figure(1)
% plot(x,H);
% axis([0 max(x) -1/c_eta_0 1/c_eta_0]);
end