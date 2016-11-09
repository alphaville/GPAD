function [x,fval,exitflag,y,output] = quadprog2(Q,q,A,bmin,bmax)

%{
n=100;m=2*n;
Q=randn(n,n);Q=Q'*Q;q=randn(n,1);A=randn(m,n);bmin=-ones(m,1);bmax=ones(m,1);
opts.Algorithm='interior-point-convex';
tic;[x,fval,exitflag,output] = quadprog(Q,q,[A;-A],[bmax;-bmin],[],[],[],[],[],opts);toc
tic;[x2,fval2,exitflag2,y2,output2] = quadprog2(Q,q,A,bmin,bmax);toc
%}

%   Copyright 2013 Panagiotis Patrinos




exitflag = 0;
fval = inf;
output = [];
m = size(A,1);
maxiter = 1e6;% maximum number of iterations
tolf    = 1e-3;% tolerance for duality gap
tolg    = 1e-3;% tolerance for feasibility
sigma   = 0.5;% reduction factor for the line search
gamma   = 100;
bmintol = bmin - tolg;
bmaxtol = bmax + tolg;
%% Offline calculations
L = chol(Q,'lower');
D = -L'\(L\A');
d = -L'\(L\q);
%% Initialize iterates
ls_iters = 0;
% initial dual iterates: y_0, y_{-1}
y  = zeros(m,1);
yp = y;
% initial primal iterates: y_0, y_{-1}
x  = d;
xp = x;
%% Ax_0, Ax_{-1}
g  = A*d;
gp = g;
% initial dual cost
dual_cost = -0.5*q'*d;
theta = 1;thetap = 1;
for k = 1:maxiter
    beta = theta*(1/thetap-1);
    w  = y + beta*(y-yp);
    xw = x + beta*(x-xp);
    gw = g + beta*(g-gp);
    % store previous values
    yp = y;
    xp = x;
    gp = g;
    thetap = theta;
    dual_costp = dual_cost;
    % calculate dual smooth cost at w
    sw  = w'*gw;
    fconjw  = -0.5*(sw+q'*xw);
    % calculate y test point
    z  = min(max(gw + w/gamma,bmin),bmax);
    y  = w + gamma*(gw-z);
    x  = D*y+d;
    g  = A*x;
    s  = y'*g;
    % calculate dual smooth cost at y
    fconj = -0.5*(s+q'*x);
    yw = y - w;
    while fconj >= fconjw - gw'*yw + 1/(2*gamma)*(yw'*yw)
        ls_iters=ls_iters+1;
        gamma = sigma*gamma;
        z  = min(max(gw + w/gamma,bmin),bmax);
        y  = w + gamma*(gw - z);
        x  = D*y + d;
        g  = A*x;
        s  = y'*g;
        yw = y - w;
        fconj  = - 0.5*(s + q'*x);
    end
    dual_cost = fconj + bmin'*min(y,0)+bmax'*max(y,0);
    fval = -s - fconj;
    gap = fval + dual_cost;
    % check for termination
    if all(g <= bmaxtol) && all(bmintol <= g)
        if  gap <= tolf
            exitflag = 1;
            if nargout>4
                output.iterations = k;
                output.constrviolation = max(max(g-bmax),max(bmin-g));
                output.gap = gap;
                output.ls_iters = ls_iters;
            end            
            return
        end
    end
    % enforce monotonicity
    if dual_cost > dual_costp
        y = yp;
        g = gp;
        x = xp;
        dual_cost = dual_costp;
    end
    theta = theta*(sqrt(theta^2+4)-theta)/2;
end
  
