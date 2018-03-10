sim = 500;
spsyc_ans=1:sim;
p = ProgressBar(sim); 
for i = 1:sim
    %pe = - (i-1)*100/(sim-1); Q4
    pe = 0.5 + (i-1)*9.5/(sim-1);
    spsyc_ans(i)=intfire(pe);
    p.progress
end
p.stop
%plot(linspace(0, -100, sim),spsyc_ans, 'LineWidth',2) Q4
plot(linspace(0.5, 20, sim),spsyc_ans, 'LineWidth',2)
xlabel('kick size (mV)') %label for 1 fig
ylabel('spikesbeforesync') %label for 1 fig

function spsyc=intfire(pe)
% simulates integrate and fire dynamics of N pulse-coupled neurons

    N = 2; % number of neurons
    
    % numerical parameters
    tf = 250; % ms
    showplot = false; %usually true, but not when running loopidy loop loop
    
    % physical parameters
    param.tau  =  10*ones(N,1); % ms, membrane time constant
    param.vr   =   0*ones(N,1); % mV, reset voltage
    param.vth  =  70*ones(N,1); % mV, threshold voltage 
    param.vl   = 100*ones(N,1); % mV, leak current reversal potential
    param.tref =   10*ones(N,1); % ms, refractory period
	
    param.epsilon = pe; % mV, kick size
    
    % connections
    % C(i,j) = 1 means a spike in neuron i kicks neuron j
    C = ones(N,N) - eye(N); % all-to-all connections *ORIGINAL VALUE%
    %C = [0,1;0,0]; %Value for Q2
    %C = zeros(N,N); C(1,2:end) = 1; % one drives the rest
    
    % initial conditions
    v0 = linspace(max(param.vr),min(param.vth),N+2)'; v0 = v0(2:N+1); % uniformly spaced between reset and threshold
    
    % plot setup
    if showplot
        clrs = jet(N);
        figure(); hold on;
        for n=1:N
            plot([0 tf] , [1 1]*param.vth(n), '--k')
            plot([0 tf] ,  [1 1]*param.vr(n), '--k')
        end
        xlabel('time (ms)')
        ylabel('membrane voltage (mV)')
        xlim([0 tf])
        ylim([min(param.vr)-5,max(param.vth)+5])
    end
    
    % events setup
    options = odeset('Event',@(t,v) event_crossthreshold(t,v,param));

    % solve the ODE
    t = 0;  % record and update the time
    v = v0; % values of the voltages
    spikesbeforesync = 0; sync = false; % record spikes until all the neurons fire together
    tref = zeros(N,1); % neurons are refractory until these times

    while t < tf
        % simulate the sub-threshold dynamics (until a neuron crosses its threshold)
        [tval,vval,te,ve,ie] = ode23(@(t,v) rhs(t,v,param,tref),[t tf],v,options);
        
        if isempty(te) % no spikes
            t = tval(end);
            spiked = false(N,1);
        else % spike
            t = te(1);                 % update the time
            v = ve(1,:)';              % update the voltages
            
            % process voltage kicks due to spikes, the spike cascade
            spiked = false(N,1);
            spiking = false(N,1);
            spiking(ie) = true;
            while any(spiking)
                n = find(spiking,1);                                 % index of a spiking neuron
                v = v + param.epsilon*C(n,:)'.*(t>tref);             % kick the neurons connected to n that are not refractory
                spiked(n) = true; spiking(n) = false;                % move n from the spiking set to the spiked set
                spiking( v>=param.vth & ~(spiking|spiked) ) = true;  % determine which neurons went over threshold
            end
            v(spiked) = param.vr(spiked);                            % reset neurons above threshold
            tref(spiked) = t + param.tref(spiked);                   % update refractory time
            
            % count spikes before sync
            sync = sync || all(spiked);
            if ~sync
                spikesbeforesync = spikesbeforesync + 1;
            end
        end
        
        % plot
        if showplot
            for n=1:N
                plot(tval,vval(:,n),'-','LineWidth',2,'Color',clrs(n,:))
                if spiked(n)
                    plot(t*[1 1],[param.vth(n),param.vr(n)],'-','LineWidth',2,'Color',clrs(n,:))
                end
            end
        end
    end
    spsyc = spikesbeforesync;
end

function dvdt = rhs(t,v,param,tref)
% tau dV/dt = Vl - V
    dvdt = (param.vl-v)./param.tau; % dynamics
    dvdt(t<tref) = 0 ; % impose refractory period
end

function [value,isterminal,direction] = event_crossthreshold(t,v,param)
    value = v - param.vth;         % watch for v = vth for each neuron
    isterminal = true(size(v));    % always stop the simulation when this occurs
    direction  = 1*ones(size(v));  % only seek up crossings
end