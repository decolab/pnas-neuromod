function [b] = BOLD(T,r)

% The Hemodynamic model with one simplified neural activity
% 
% BOLD(f0,T,n_fig,suite)
%
% T       : total time (s)

%global itaus itauf itauo ialpha Eo dt

ch_int = 0;         % 0: Euler, 1: ode45

dt  = 0.001;        % (s)
t0  = (0:dt:T)';
n_t = length(t0);

t_min = 20;
n_min = round(t_min/dt);

r_max = max(r);

% BOLD model parameters

taus   = 0.65; % 0.8;    % time unit (s)
tauf   = 0.41; % 0.4;    % time unit (s)
tauo   = 0.98; % 1;      % mean transit time (s)
alpha  = 0.32; % 0.2;    % a stiffness exponent
itaus  = 1/taus;
itauf  = 1/tauf;
itauo  = 1/tauo;
ialpha = 1/alpha;
Eo     = 0.34; % 0.8;    % resting oxygen extraction fraction
vo     = 0.02;
k1     = 7*Eo; 
k2     = 2; 
k3     = 2*Eo-0.2;

Eo=0.4;
TE=0.04;
vo=0.04;
k1=4.3*40.3*Eo*TE;
k2=25*Eo*TE;
k3=1;

% Initial conditions

x0  = [0 1 1 1];

tic;

if ch_int == 0
    
    % Euler method

    t      = t0;
    x      = zeros(n_t,4);
    x(1,:) = x0;
    for n = 1:n_t-1;
        x(n+1,1) = x(n,1) + dt*( r(n)-itaus*x(n,1)-itauf*(x(n,2)-1) );
        x(n+1,2) = x(n,2) + dt*x(n,1);
        x(n+1,3) = x(n,3) + dt*itauo*(x(n,2)-x(n,3)^ialpha);
        x(n+1,4) = x(n,4) + dt*itauo*(x(n,2)*(1-(1-Eo)^(1/x(n,2)))/Eo - (x(n,3)^ialpha)*x(n,4)/x(n,3));
    end

else
        
    %opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
    opt = odeset('RelTol',1e-6,'AbsTol',1e-6);

    [t,x] = ode45('BOLD_ODEs',t0,x0,opt);

end

t  = t(n_min:end);
s  = x(n_min:end,1);
fi = x(n_min:end,2);
v  = x(n_min:end,3);
q  = x(n_min:end,4);
%b  = 100/Eo*vo*( k1.*(1-q) + k2*(1-q./v) + k3*(1-v) );
b  = vo*( k1.*(1-q) + k2*(1-q./v) + k3*(1-v) );
clear x;





