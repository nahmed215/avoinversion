function [ x, g, f, info, count, img] = linesearch_more(fg, xp, f, g, s, stp, ftol,gtol, max_linesearch, min_step, max_step)
x=xp;
count=0;

%     /* Check the input parameters for errors. */
if (stp <= 0)
    info=999;
    return;
end

%     /* Compute the initial gradient in the search direction. */
dginit=sum(g.*s);

%     /* Make sure that s points to a descent direction. */
if (dginit > 0)
    info=998;
    return;
end

% Initialize local variables. */
brackt = 0;
stage1 = 1;
finit = f;
dgtest = ftol * dginit;
width = max_step - min_step;
prev_width = 2.0 * width;

%
%     /*
%         The variables stx, fx, dgx contain the values of the step,
%         function, and directional derivative at the best step.
%         The variables sty, fy, dgy contain the value of the step,
%         function, and derivative at the other endpoint of
%         the interval of uncertainty.
%         The variables stp, f, dg contain the values of the step,
%         function, and derivative at the current step.
%     */
stx = 0;
sty = 0.;
fx = finit;
fy = finit;
dgx = dginit;
dgy = dginit;
infou=0;
while (1)
    
    %         /*
    %             Set the minimum and maximum steps to correspond to the
    %             present interval of uncertainty.
    %          */
    if (brackt)
        stmin = min(stx, sty);
        stmax = max(stx, sty);
    else
        stmin = stx;
        stmax = stp + 4.0 * (stp - stx);
    end
    
    %         /* Clip the step in the range of [stpmin, stpmax]. */
    if (stp < min_step)
        stp = min_step;
    end
    if (max_step < stp)
        stp = max_step;
    end
    
    %         /*
    %             If an unusual termination is to occur then let
    %             stp be the lowest point obtained so far.
    %          */
    if ((brackt && ((stp <= stmin || stmax <= stp) || max_linesearch <= count + 1 || infou ~= 0)) || (brackt && (stmax - stmin <= eps * stmax)))
        stp = stx;
    end
    
    %         /*
    %             Compute the current value of x:
    %                 x <- x + (*stp) * s.
    %          */
    x=xp;
    x=x+ stp*s;
    
    %         /* Evaluate the function and gradient values. */
    if (nargout>5)
        [f,g,img]=feval(fg,x);
    else
        [f,g]=feval(fg,x);
    end
    dg=sum(g.*s);
    
    
    
    ftest1 = finit + stp * dgtest;
    count=count+1;
    
    %         /* Test for errors and convergence. */
    if (brackt && ((stp <= stmin || stmax <= stp) || infou ~= 0))
        %             /* Rounding errors prevent further progress. */
        info=997;
        return;
    end
    if (stp == max_step && f <= ftest1 && dg <= dgtest)
        %              /* The step is the maximum value. */
        info=996;
        return;
    end
    if (stp == min_step && (ftest1 < f || dgtest <= dg))
        %             /* The step is the minimum value. */
        info=995;
        return;
    end
    if (brackt && (stmax - stmin) <= eps * stmax)
        %             /* Relative width of the interval of uncertainty is at most xtol. */
        info=994;
        return;
    end
    if (max_linesearch <= count)
        %             /* Maximum number of iteration. */
        info=993;
        return;
    end
    if (f <= ftest1 && abs(dg) <= gtol *(-dginit))
        %             /* The sufficient decrease condition and the directional derivative condition hold. */
        info=0;
        return;
    end
    
    %         /*
    %             In the first stage we seek a step for which the modified
    %             function has a nonpositive value and nonnegative derivative.
    %          */
    if (stage1 && f <= ftest1 && min(ftol, gtol) * dginit <= dg)
        stage1 = 0;
    end
    
    %         /*
    %             A modified function is used to predict the step only if
    %             we have not obtained a step for which the modified
    %             function has a nonpositive function value and nonnegative
    %             derivative, and if a lower function value has been
    %             obtained but the decrease is not sufficient.
    %          */
    if (stage1 && ftest1 < f && f <= fx)
        %             /* Define the modified function and derivative values. */
        fm = f - stp * dgtest;
        fxm = fx - stx * dgtest;
        fym = fy - sty * dgtest;
        dgm = dg - dgtest;
        dgxm = dgx - dgtest;
        dgym = dgy - dgtest;
        
        %             /*
        %                 Call update_trial_interval() to update the interval of
        %                 uncertainty and to compute the new step.
        %              */
        [stx, fxm, dgxm, sty, fym, dgym, stp, brackt, infou] = update_trial_int(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, stmin, stmax, brackt);
        
        %             /* Reset the function and gradient values for f. */
        fx = fxm + stx * dgtest;
        fy = fym + sty * dgtest;
        dgx = dgxm + dgtest;
        dgy = dgym + dgtest;
    else
        %             /*
        %                 Call update_trial_interval() to update the interval of
        %                 uncertainty and to compute the new step.
        %              */
        [stx, fx, dgx,sty, fy, dgy, stp, brackt, infou] = update_trial_int(stx, fx, dgx,sty, fy, dgy,stp, f, dg,stmin, stmax, brackt);
    end
    
    %         /*
    %             Force a sufficient decrease in the interval of uncertainty.
    %          */
    if (brackt)
        if (0.66 * prev_width <= abs(sty - stx))
            stp = stx + 0.5 * (sty - stx);
        end
        prev_width = width;
        width = abs(sty - stx);
    end
end
info=993;
return;
end



