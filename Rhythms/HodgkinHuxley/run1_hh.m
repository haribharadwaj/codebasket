
%% Script to 'burn in' an initial condition
% Problem 1B
Tsolve = [0 50];
init_cond = [0; 0; 1; 0.5];
Istim = 0;

[t,p] = ode45(@(t,p) hh_rhs(t,p,Istim),Tsolve,init_cond);

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
init_cond = [-8.4905e-04; 0.0529; 0.5975; 0.3177];

%% Real simulations

% Problem 1C
Tsolve = [0 50];


% Determining threshold (for 50ms long simulation)
Ilist = 0:-0.01:-2; % 2 decimals accuracy, -ve is stimulating => -0.01
for j = 1:numel(Ilist)
    [t,p] = ode45(@(t,p) hh_rhs(t,p,Ilist(j)),Tsolve,init_cond);
    if(any(p(:,1) < -80))
       Ithresh_50ms = Ilist(j);
       break;
    end
end

disp(Ithresh_50ms);
[t,p] = ode45(@(t,p) hh_rhs(t,p,Ithresh_50ms),Tsolve,init_cond);
figure;
plot(t,p(:,1));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Membrane Potential (mV)','FontSize',14,'FontWeight','bold');

% The threshold was -1.95 microA

% Subthreshold (just) stimulation
Istim = Ithresh_50ms + 0.01;
[t,p] = ode45(@(t,p) hh_rhs(t,p,Istim),Tsolve,init_cond);
figure;
plot(t,p(:,1));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Membrane Potential (mV)','FontSize',14,'FontWeight','bold');


% The depolarization for 0.01 mV 'short of' the threshold  is very small
% (about 8 mV) suggesting that it is a true 'all or nothing' response
% suggesting that it is an annihilation of a fixed point pair. Thus it
% suggests that it is a SADDLE NODE BIFURCATION on a cycle. If it were a
% hopf kind of route to a cycle, we would expect the amplitude of cycle to
% grow slowly with change is parameter value.


%% Problem 1D
figure;
for Istim = -1.95:-20:-400
    [t,p] = ode45(@(t,p) hh_rhs(t,p,Istim),Tsolve,init_cond);
    plot(t,p(:,1));
pause;
end
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
ylabel('Membrane Potential (mV)','FontSize',14,'FontWeight','bold');


% Looking at the plots for a range from the threshold of about -2 upto -400
% it seems like the frequency increases continuously with Istim
% (magnitude) and the amplitude of the action potential seems to decrease
% continously ang get damped out gradually to a flat response. I could not
% identify any particular threshold above which this happens.. It seems
% like the after a point, the amplitude just gradually drops to zero.


%% Problem 3
% The real dimensionality of the model...

Istim = -3; % Arbitrary suprathreshold stimulus current

[t,p] = ode45(@(t,p) hh_rhs(t,p,Istim),Tsolve,init_cond);
figure;
subplot(2,3,1);
plot(p(:,1),p(:,2));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('V(t) (ms)','FontSize',14,'FontWeight','bold');
ylabel('m(t)','FontSize',14,'FontWeight','bold');

subplot(2,3,2);
plot(p(:,1),p(:,3));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('V(t) (ms)','FontSize',14,'FontWeight','bold');
ylabel('h(t)','FontSize',14,'FontWeight','bold');

subplot(2,3,3);
plot(p(:,1),p(:,4));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('V(t) (ms)','FontSize',14,'FontWeight','bold');
ylabel('n(t)','FontSize',14,'FontWeight','bold');

subplot(2,3,4);
plot(p(:,2),p(:,3));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('m(t)','FontSize',14,'FontWeight','bold');
ylabel('h(t)','FontSize',14,'FontWeight','bold');

subplot(2,3,5);
plot(p(:,2),p(:,4));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('m(t)','FontSize',14,'FontWeight','bold');
ylabel('n(t)','FontSize',14,'FontWeight','bold');

subplot(2,3,6);
plot(p(:,3),p(:,4));
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('h(t)','FontSize',14,'FontWeight','bold');
ylabel('n(t)','FontSize',14,'FontWeight','bold');


% It seems like on the limit cycle, the value of n(t) and h(t) are
% approximately related by a linear curve with h-intercept 0.6 and
% n-intercept about 0.8: h(t) = -(n(t) - 0.8)*(0.6/0.8) => h =
% (0.8-n)*0.75. This has been used for problem 3B

%% 3D simplification (Problem 3B)
%h is eliminates
Tsolve = [0 50];
init_cond_3D = [-8.4905e-04; 0.0529; 0.3177];

% burn in
Istim = 0;
[t,p3D] = ode45(@(t,p3D) hh_rhs_3D(t,p3D,Istim),Tsolve,init_cond_3D);

init_cond_3D = [p3D(end,1); p3D(end,2); p3D(end,3)];
Istim = -10; %microV
[t,p3D] = ode45(@(t,p3D) hh_rhs_3D(t,p3D,Istim),Tsolve,init_cond_3D);
figure;
plot(t,p3D(:,1));
set(gca,'FontSize',14,'FontWeight','bold');
ylabel('Membrane Voltage (mV)','FontSize',14,'FontWeight','bold');
xlabel('Time (ms)','FontSize',14,'FontWeight','bold');
% One is able to get the same response albeit at stronger (more negative)
% stimulus currents. It seems to capture the overall shape of the action
% potential and subthreshold trajectories reasonably well. The resting
% potential for the unforced case (Istim = 0) seems to be slightly altered
% as well.


%% Covariation of threshold with amplitude and duration of Istim pulse
% Problem 4
% A simple squarePulse.m function has been defined to specify Istim as a
% function of time
Tsolve = [0 500];
duration = 1:10:500;
Ithresh = zeros(size(duration));
for j = 1:numel(duration);
    for Imax = 0:-0.2:-4
        [t,p] = ode45(@(t,p) hh_rhs(t,p,squarePulse(t,Imax,duration(j))),Tsolve,init_cond);
        if(any(p(:,1) < -80))
            Ithresh(j) = Imax;
            break;
        end
    end
end

plot(duration,Ithresh);
set(gca,'FontSize',14,'FontWeight','bold');
title('Covariation of threshold with amplitude and duration of Istim pulse','FontSize',14,'FontWeight','bold');
ylabel('Threshold Current for Given Duration (microA)','FontSize',14,'FontWeight','bold');
xlabel('Pulse Duration (ms)','FontSize',14,'FontWeight','bold');


% It seems like if one waits long enough, even with a weaker pulse, one is
% able to get an action potential and vice versa.


        
        
