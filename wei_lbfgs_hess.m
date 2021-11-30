function [xout,pf,iter,y,s,YS,alpha,fim,bound]=wei_lbfgs_hess(fg, x, m, max_iterations, linesearch_type, max_linesearch, ftol, gtol, wolfe, delta, pgtol, fhess)

%     Constant parameters and their default values. 
if (nargin <3)
    m=9;
end
if (nargin <4)
    max_iterations=0;
end
if (nargin <5)
    linesearch_type=4;
end
if (nargin <6)
    max_linesearch=10;
end
if (nargin <7)
    ftol=1e-4;
end
if (nargin <8)
    gtol=0.9;
end
if (nargin <9)
    wolfe=1e-3;
end
if (nargin <10)
    delta=1e-2;
end
if (nargin <11)
    pgtol=1e-5;
end


n=length(x);
s=zeros(n,m);
y=zeros(n,m);
alpha=zeros(m,1);
xout=zeros(n, 1);
%     Evaluate the function value and its gradient. 
[fx,g]=feval(fg,x);
nfg=1;
%     Store the initial value of the objective function. 

pf(1) = fx;


%         Compute the direction we assume the initial hessian matrix H_0 as
%         the identity matrix
d=-g;

%        Make sure that the initial variables are not a minimizer.

xnorm=sqrt(sum(x.*x));
gnorm=sqrt(sum(g.*g));

if (xnorm < 1.0)
    xnorm = 1.0;
end
if (gnorm / xnorm <= pgtol)
    display('LBFGS_ALREADY_MINIMIZED');
    iter=1;
    xout=x;
    return;
end

%     Compute the initial step:
step = 1.0 / sqrt(sum(d.* d));

steps(1) = step;

iter = 1;
fim = 0;
while (1)
%     clc;
    save pre_linesearch.mat x g;
    display(iter);
    display(nfg);
    %         Store the current position and gradient vectors. 
    xp=x;
    gp=g;
    
    %          Search for an optimal step. 
    if(linesearch_type<4)
        [x, g, fx, ls, count]=linesearch_backtrack(fg, x,  fx, g, d, step, ftol,wolfe, max_linesearch, 1e-20, 1e20,linesearch_type);
    else
        [x, g, fx, ls, count]=linesearch_more(fg, x, fx, g, d, step, ftol, gtol, max_linesearch, 1e-20, 1e20);
    end
    nfg=nfg+count;
    

    if (ls ~= 0)
%       Revert to the previous point.
        iter=iter-1;
        x=xp;
        g=gp;
        
        display('Could not find accetaple step, finishing.');                
        return
    end    
    %        Compute x and g norms. 
    xnorm=sqrt(sum(x.*x));
    gnorm=sqrt(sum(g.*g));
    
    
    % Report the progress.
%              Store the current value of the objective function. 
    pf(iter)= fx;
    steps(iter) = step;
    if(iter==1)
        xout=x;
    else
        xout=[xout, x];
    end
    save pre_results.mat xout nfg;
    
    %             Convergence test.
    %             The criterion is given by the following formula:
    %                 |g(x)| / \max(1, |x|) < \epsilon
    
    if (xnorm < 1.0)
        xnorm = 1.0;
    end
    if (gnorm / xnorm <= pgtol)
        display('Function MINIMIZED by gnorm/xnorm');        
        return;
    end
     %             Convergence test.
    %             The criterion is given by the following formula:
    %             (fx(i+1)-fx(i))/fx(1)
    if(iter>1)
        dfx=abs((pf(iter-1)-fx)/fx);
            if (dfx<delta)
                display('Function MINIMIZED by difference in objective function');
                return;
            end
    end
    
    

    
    if (max_iterations ~= 0 && max_iterations < iter+1)
        % Maximum number of iterations.
        display('Maximum number of iterations');
        return;
    end
    
    
    %             Update vectors s and y:
    %                 s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
    %                 y_{k+1} = g_{k+1} - g_{k}.
    
    s(:,fim+1)=x-xp;
    y(:,fim+1)=g-gp;
    
    
    %             Compute scalars ys and yy:
    %                 ys = y^t \cdot s = 1 / \rho.
    %                 yy = y^t \cdot y.
    %             Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
    
    ys=sum(s(:,fim+1).*y(:,fim+1));
    yy=sum(y(:,fim+1).*y(:,fim+1));
    YS(fim+1)=ys;
    
    
    %        
    %             Recursive formula to compute dir = -(H \cdot g).
    %                 This is described in page 779 of:
    %                 Jorge Nocedal.
    %                 Updating Quasi-Newton Matrices with Limited Storage.
    %                 Mathematics of Computation, Vol. 35, No. 151,
    %                 pp. 773--782, 1980.
    %          
    if(m <= iter)
        bound= m;
    else
        bound= iter;
    end
    iter=iter+1;
    fim = mod((fim + 1), m);
    
    % Compute the descent direction.
    
    % Compute the negative of gradients.
    d=-g;
    
    j = fim;
    for i = 1:bound
        j = mod((j + m - 1),m);  %  /* if (--j == -1) j = m-1;
        %\alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}.
        alpha(j+1)=sum(s(:,j+1).*d);
        alpha(j+1)= alpha(j+1)/YS(j+1);
        %q_{i} = q_{i+1} - \alpha_{i} y_{i}. */
        d=d+ (y(:,j+1)*-alpha(j+1));
    end
    d=d*ys/yy;
           
    for i = 1:bound
        %\beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}.
        beta=sum(y(:,j+1).*d);
        beta =beta/ YS(j+1);
        %\gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}. 
        d=d+s(:,j+1)*(alpha(j+1)-beta);
        j = mod(j + 1, m);    %     if (++j == m) j = 0; 
    end
    
% Now the search direction d is ready. We try step = 1 first.
  
    step = 1.0;
end

% %  Return the final value of the objective function. *

% pf(iter+1)=fx;
% display('finished')
end