function f = force(t,Amplitude,omega,tMax);
if t < tMax/2
    f = Amplitude;
else 
%     f = Amplitude*t;
    f = Amplitude*cos(omega*(t-tMax/2));
end
end