function dxdt = HHrhs(t,x)

    setparameters;
    setgatevars;
    
    v = x(1);
    q = x(2:4)';
    
    % original model (a)
    g =  g_Na*q(mi)^3*q(hi)      + g_K*q(ni)^4     + g_L;
    Eg = g_Na*q(mi)^3*q(hi)*E_Na + g_K*q(ni)^4*E_K + g_L*E_L;
    
    % you need to modify the expressions above for each approximation
    % for example to impose m=m_inf, use qinf{mi}(v) instead of q(mi)
    
    dxdt = NaN(4,1);
    dxdt(1) = -g*v/C + Eg/C;
    
    for gv=1:ngv
        dxdt(1+gv) = alpha{gv}(v)*(1-q(gv)) - beta{gv}(v)*q(gv);
    end
end