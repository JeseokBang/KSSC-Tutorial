function P_Coeff=AssoLegen_Coeff(c_val,m,l) %P_(l)^(m)
s_val=sqrt(1-c_val^2);
P_Coeff=0;
for k=m:l
    C1=gamma(k+1)/gamma(k-m+1);
    C2=gamma(l+1)/gamma(k+1)/gamma(l-k+1)*gamma((l+k-1)/2+1)/gamma(l+1)/gamma((l+k-1)/2-l+1);
    C3=C1*C2;
    P_Coeff=P_Coeff+C3*c_val^(k-m);
end
P_Coeff=(-1)^m*2^l*s_val^m*(P_Coeff);
end