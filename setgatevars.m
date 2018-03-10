% defines the alpha, beta functions for the gating variables
% the function handles are stored in cell arrays
% alpha{1}(V) computes alpha_m of each value in the array V
ngv = 3; % number of gating variables
mi=1;ni=2;hi=3;
alpha = cell(ngv,1);
alpha{mi} = @(V) (0.1*V+4.5)./(1-exp(-0.1*V-4.5));
alpha{ni} = @(V) 0.1*(0.1*V+6)./(1-exp(-0.1*V-6));
alpha{hi} = @(V) 0.07*exp(-0.05*V-3.5);

beta = cell(ngv,1);
beta{mi} = @(V) 4*exp(-(V+70)/18);
beta{ni} = @(V) 0.125*exp(-(V+70)/80); 
beta{hi} = @(V) 1./(1+exp(-0.1*V-4));

qinf = cell(ngv,1);
for gv=1:ngv
    qinf{gv} = @(V) alpha{gv}(V)./(alpha{gv}(V)+beta{gv}(V));
end