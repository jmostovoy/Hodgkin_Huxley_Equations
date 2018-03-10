% solves (numerically) the cable equation
clear;

% setup parameters and gating functions
setparameters;
setgatevars;

% numerical parameters
N        = 2^10;       %          number of grid points
klokmax  = 6500;       %          number of time steps
tf       = 65;         % ms       total simulation time
dt       = tf/klokmax; % ms       timestep

% display parameters
lw = 2;                      % line width
sty = {'-r','-b','-k'};      % gating variable styles
frameevery = 10;             % show a frame every frameevery timesteps

% setup
dx = L/N;
x = (0:N-1)'*dx;

% initial condition
% solve for the resting membrane potential
Vr = fzero(@(vr) (E_Na*g_Na*qinf{mi}(vr)^3*qinf{hi}(vr) + E_K*g_K*qinf{ni}(vr)^4+E_L*g_L) / (g_Na*qinf{mi}(vr)^3*qinf{hi}(vr) + g_K*qinf{ni}(vr)^4+g_L) - vr,-70);
q = ones(N,ngv);
for gv=1:ngv
    q(:,gv) = qinf{gv}(Vr);
end

%V = Vr*ones(N,1);
V = Vr + abs(Vr)*exp(-(x-L/5).^2/0.5^2); % increase the voltage near 1/5 the length

% setup for solving the cyclic tridiagonal system
% only the main triagonal changes from timestep to timestep
sd  = -a/(4*rho)*dt/dx^2/C*ones(N,1); % sub/super diagonal value
A = sparse(N,N);
A(2:N+1:end) = sd(2:N);
A(N+1:N+1:end) = sd(1:N-1);
A(1,N) = sd(1);
A(N,1) = sd(N);

% plotting setup
fh = figure();
frame = 0;

% setup for estimating the wave speed
maxV = -Inf;
maxt = NaN;
ind = round(0.5*N);

% step through time
for klok=1:klokmax
    
    % update the gating variables
    % (gating variables are defined at half timesteps)
    for gv=1:ngv
        alphavals = alpha{gv}(V);
        betavals  =  beta{gv}(V);
        q(:,gv)   = ( q(:,gv).*( 1 - dt/2*(alphavals+betavals) ) + dt*alphavals ) ./ (1 + dt/2*(alphavals+betavals));
    end
    
    % update the potential
    % by solving a cyclic tridiagonal system
    g =  g_Na*q(:,mi).^3.*q(:,hi)      + g_K*q(:,ni).^4     + g_L;
    Eg = g_Na*q(:,mi).^3.*q(:,hi)*E_Na + g_K*q(:,ni).^4*E_K + g_L*E_L;
    
    A(1:N+1:end)  = 1 + g/C*dt/2 + a/(2*rho)*dt/dx^2/C*ones(N,1); % main diagonal value
    
    rhs =   V .* (1 - g/C*dt/2 - a/(2*rho)*dt/dx^2/C ) ...
          + a/(4*rho)*dt/dx^2* ( V([2:N,1]) + V([N,1:(N-1)]) )/C ...
          + Eg/C*dt; % right hand side

    V = A\rhs;
    
    % plot
    if mod(klok,frameevery) == 0
        
        frame = frame + 1;
        figure(fh); clf;
    
        subplot(2,1,1); hold on;
        plot([x;L],[V;V(1)],'-k','LineWidth',lw)
        axis([0 L -100 40])
        title(['time=' num2str(klok*dt,'%10.3f') ' ms'])
        grid on
        ylabel('V (mV)')
        plot(x(ind),maxV,'or')
        
        subplot(2,1,2); hold on;
        for gv=1:ngv
            plot([x;L],[q(:,gv);q(1,gv)],sty{gv},'LineWidth',lw)
        end
        axis([0 L 0 1])
        grid on
        xlabel('x (cm)')
        ylabel('m, n, h')

        
        drawnow;
        
    end

    % estimate the AP speed
    if V(ind)>maxV
        maxV = V(ind);
        maxt = klok*dt;
    end
    
end % of timestepping

display([ 'AP speed: ' num2str( (L/2-L/5)/maxt ) ' cm/ms']);