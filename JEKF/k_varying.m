function k_t = k_varying(k_C, k_P, omega,t,tMax)

tChange = tMax/10;
if (t < tChange)
    k_t = k_C;
else
%     k_t = k_C + k_P;
    k_t = k_C + k_P*sin(omega*(t-tChange));
end

end