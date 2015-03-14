function dpdt = hh_rhs(t,p,Istim) %#ok<INUSL>

% Function that implements the RHS of the Hodgkin-Huxley model. All
% parameters are hard-coded. Istim is made an additional parameter.
%
% Usage:
% dpdt = hh_rhs(t,p,Istim);
%   t: The time instant (scalar)
%   p: The set of state variables: V,m,h,n (In that order in a column)
%   Istim: Input stimulus current
%
% To use with ODE45 create an anonymous function with just 2 variables t,p
% from hh_rhs() just like you would do with some inbuilt matlab function.
%
% Example: [t,p] = ode45(@(t,p) hh_rhs(t,p,-2),[0 50],[0;0;1;0.5]);
%   This is similar to [t,x] = ode45(@(t,x) normpdf(x,0,1),[0 5],1);
% One could even make Istim a function of t.
%------------
% Hari Bharadwaj, Dec 4, 2010
%------------


V = p(1);
m = p(2);
h = p(3);
n = p(4);

% Maximal conductances in mmho/cm^2
GNa = 120;
GK = 36;
Gl = 0.3;

% Reversal potentials for the different ion species in mV
ENa = -115;
EK = 12;
El = -10.5989;

% Membrance capacitance microF/cm^2
C = 0.775;


alpha_m = 0.1*(V+25)/(exp((V+25)/10) - 1);
beta_m = 4*exp(V/18);
alpha_h = 0.07*exp(V/20);
beta_h = 1/(exp((V+30)/10) + 1);
alpha_n = 0.01*(V+10)/(exp((V+10)/10) - 1);
beta_n = 0.125*exp(V/80);


gNa = GNa*m^3*h;
gK = GK*n^4;
gl = Gl;

INa = gNa*(ENa - V);
IK = gK*(EK - V);
Il = gl*(El - V);

% The dynamic terms:

dVdt = (INa + IK + Il + Istim)/C; % Note: Istim is an input paramter
% The units above are self consistent mV/ms
% I would hence be in microA

dmdt = alpha_m*(1-m) - beta_m*m;
dhdt = alpha_h*(1-h) - beta_h*h;
dndt = alpha_n*(1-n) - beta_n*n;

dpdt = [dVdt; dmdt; dhdt; dndt];

