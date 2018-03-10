% solves the space-clamped Hodgkin-Huxley equations
%clear; Was here originally, was not a fan
sims = 100;
V0 = 1:sims;
v = ones(16384, sims);
p = ProgressBar(sims); 
for i  = 1:sims

% setup parameters and gating functions
setparameters;
setgatevars;

tf = 50;    % ms    final time
npts = 2^14; %       number of points in the computed solution

% solve for the resting membrane potential
Vr = fzero(@(vr) (E_Na*g_Na*qinf{mi}(vr)^3*qinf{hi}(vr) + E_K*g_K*qinf{ni}(vr)^4+E_L*g_L) / (g_Na*qinf{mi}(vr)^3*qinf{hi}(vr) + g_K*qinf{ni}(vr)^4+g_L) - vr,-70);

% setup the initial condition
q0 = ones(1,ngv);
for gv=1:ngv
    q0(gv) = qinf{gv}(Vr); % set each gating variable to its value at the resting membrane potential
end
V0(i) = -sims/2 + (i-1)*sims/(sims-1); % mV

% solve the ode
[t,vars] = ode45( @HHrhs, linspace(0,tf,npts), [V0(i),q0] );

% extract the solution
v(:,i) = vars(:,1);
m = vars(:,mi+1);
n = vars(:,ni+1);
h = vars(:,hi+1);

% plot
wannaplot = 0;
if wannaplot == 1
    figure(); hold off;
    plot(t,v,'-k','LineWidth',2)
    xlabel('t (ms)')
    ylabel('V (mV)')

    % TODO: plot the gating varialbes m, n, and h vs t
    plot(t, [m,n,h], 'LineWidth',3)
    xlabel('Time') %label for 1 fig
    ylabel('Gating Variable') %label for 1 fig
    legend('m(t)','n(t)', 'h(t)'...
        ,'Location','northeast')
end
p.progress
end
p.stop
ColorSet = varycolor(sims);
hold on
set(gca, 'ColorOrder', ColorSet);
for k = 1:sims
plot(t, v(:,k), 'LineWidth',2, 'Color', ColorSet(k,:))
end
xlabel('Time') %label for 1 fig
ylabel('V') %label for 1 fig
lgnd_str = repmat("hi", [sims, 1]);
for j = 1:sims
    lgnd_str(j) = strcat('V0= ', num2str(V0(j)));
end
%legend(lgnd_str,'Location','northeast')
legend off
set(gcf, 'Colormap', ColorSet);
colorbar('Ticks',linspace(0,1, sims),...
         'TickLabels',char(lgnd_str))