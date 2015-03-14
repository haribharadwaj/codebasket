
%% Script to test the simplified 3D model in Problem 3B


Tsolve = [0 50];
init_cond = [0; 0; 0.5];
Istim = 4; % Hyperpolarization current

[t,p] = ode45(@(t,p) hh_rhs_3D(t,p,Istim),Tsolve,init_cond);

% Note that the call is slightly different because of that way I
% parametrized hh_rhs to include Istim as an input.

% Plotting
figure;
subplot(2,1,1)
plot(t,p(:,1)); %V(t)
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Membrane Potential (mV)','FontSize',14,'FontWeight','bold');
subplot(2,1,2)
plot(t,p(:,2),t,p(:,3),t,p(:,4)); %m,h,n as function of t
set(gca,'FontSize',14,'FontWeight','bold');
legend('m','h','n');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Gating Variables','FontSize',14,'FontWeight','bold');

%To view and manually enter the burn in for further simulations
disp(p(end,1));
disp(p(end,2:4));

% The 'RESTING' membrane voltage is now close to 5mV as opposed to zero
% 4.

%% Real simulations

% CASE 1: Resting condition burnt in with Istim = 4 micro A
Tsolve = [0 50];
init_cond = [4.9266; 0.0292; 0.7532; 0.2456];
Istim = 0;
[t,p] = ode45(@(t,p) hh_rhs(t,p,Istim),Tsolve,init_cond);
figure;
subplot(2,1,1)
plot(t,p(:,1));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Membrane Potential (mV)','FontSize',14,'FontWeight','bold');
subplot(2,1,2)
plot(t,p(:,2),t,p(:,3),t,p(:,4)); %m,h,n as function of t
set(gca,'FontSize',14,'FontWeight','bold');
legend('m','h','n');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Gating Variables','FontSize',14,'FontWeight','bold');

% CASE 2: 
% Old init_cond for m,h,n, but new V(0)
Tsolve = [0 50];
init_cond_old_mhn = [4.9266; 0.0529; 0.5975; 0.3177];
[t,p] = ode45(@(t,p) hh_rhs(t,p,Istim),Tsolve,init_cond_old_mhn);
figure;
subplot(2,1,1)
plot(t,p(:,1));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Membrane Potential (mV)','FontSize',14,'FontWeight','bold');
subplot(2,1,2)
plot(t,p(:,2),t,p(:,3),t,p(:,4)); %m,h,n as function of t
set(gca,'FontSize',14,'FontWeight','bold');
legend('m','h','n');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Gating Variables','FontSize',14,'FontWeight','bold');


% There is an action potential for CASE 1 where the stimulation was an
% increase in (negative) current of 4 microA that stimulates the neuron to
% go away from the fixed point where it is because the CASE is a natural
% rest state that has been burnt in. In CASE 2, the initial condition given
% is not close to a fixed point for that particular value of the current
% parameter. The inactivation variable 'h' is still high and the the
% initial condition is far from the separatrix that divides the two domains
% of attraction (FP and limit cycle) analogous to the Fitzhugh-Nagumo
% model. Physiologically, the cell is hyperpolarized in case 1 and the
% current depolarizes it. In case 2 however, the excess charge accumulation
% is not balanced by the chemical gradient to start with.


