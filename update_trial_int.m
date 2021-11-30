function [x, fx, dx, y, fy, dy, t, brackt, info ] = update_trial_int(x, fx, dx, y, fy, dy, t, ft, dt, tmin, tmax, brackt)


%     /* Check the input parameters for errors. */
dsign=((dt*dx/abs(dx))<0);
if (brackt)
    if (t <= min(x, y) || max(x, y) <= t)
        %             /* The trival value t is out of the interval. */
        info=999;
        return;
    end
    if (0. <= dx * (t - x))
        %             /* The function must decrease from x. */
        info=998;
        return;
    end
    if (tmax < tmin)
        %             /* Incorrect tmin and tmax specified. */
        info=997;
        return;
    end
end

done=0;
%     /*
%         Trial value selection.
%      */
if (fx < ft)
    %         /*
    %             Case 1: a higher function value.
    %             The minimum is brackt. If the cubic minimizer is closer
    %             to x than the quadratic one, the cubic one is taken, else
    %             the average of the minimizers is taken.
    %          */
    brackt = 1;
    bound = 1;
    mc=cubic_min1(x, fx, dx, t, ft, dt);
    mq=quad_min1(x, fx, dx, t, ft);
    if (abs(mc - x) < abs(mq - x))
        newt = mc;
    else
        newt = mc + 0.5 * (mq - mc);
    end
    done=1;
end
if (dsign && ~done)
    %         /*
    %             Case 2: a lower function value and derivatives of
    %             opposite sign. The minimum is brackt. If the cubic
    %             minimizer is closer to x than the quadratic (secant) one,
    %             the cubic one is taken, else the quadratic one is taken.
    %          */
    brackt = 1;
    bound = 0;
    mc=cubic_min1(x, fx, dx, t, ft, dt);
    mq=quad_min2(x, dx, t, dt);
    if (abs(mc - t) > abs(mq - t))
        newt = mc;
    else
        newt = mq;
    end
    done=1;
end
if (abs(dt) < abs(dx) && ~done)
    %         /*
    %             Case 3: a lower function value, derivatives of the
    %             same sign, and the magnitude of the derivative decreases.
    %             The cubic minimizer is only used if the cubic tends to
    %             infinity in the direction of the minimizer or if the minimum
    %             of the cubic is beyond t. Otherwise the cubic minimizer is
    %             defined to be either tmin or tmax. The quadratic (secant)
    %             minimizer is also computed and if the minimum is brackt
    %             then the the minimizer closest to x is taken, else the one
    %             farthest away is taken.
    %          */
    bound = 1;
    mc=cubic_min2(x, fx, dx, t, ft, dt, tmin, tmax);
    mq=quad_min2(x, dx, t, dt);
    if (brackt)
        if (abs(t - mc) < abs(t - mq))
            newt = mc;
        else
            newt = mq;
        end
    else
        if (abs(t - mc) > abs(t - mq))
            newt = mc;
        else
            newt = mq;
        end
    end
    done=1;
end

if (~done)
    %         /*
    %             Case 4: a lower function value, derivatives of the
    %             same sign, and the magnitude of the derivative does
    %             not decrease. If the minimum is not brackt, the step
    %             is either tmin or tmax, else the cubic minimizer is taken.
    %          */
    bound = 0;
    if (brackt)
        newt=cubic_min1(t, ft, dt, y, fy, dy);
    else
        if (x < t)
            newt = tmax;
        else
            newt = tmin;
        end
    end
end

%     /*
%         Update the interval of uncertainty. This update does not
%         depend on the new step or the case analysis above.
%
%         - Case a: if f(x) < f(t),
%             x <- x, y <- t.
%         - Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
%             x <- t, y <- y.
%         - Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0,
%             x <- t, y <- x.
%      */
if (fx < ft)
    %         /* Case a */
    y = t;
    fy = ft;
    dy = dt;
else
    %         /* Case c */
    if (dsign)
        y = x;
        fy = fx;
        dy = dx;
    end
    %         /* Cases b and c */
    x = t;
    fx = ft;
    dx = dt;
end

%     /* Clip the new trial value in [tmin, tmax]. */
if (tmax < newt)
    newt = tmax;
end
if (newt < tmin)
    newt = tmin;
end

%     /*
%         Redefine the new trial value if it is close to the upper bound
%         of the interval.
%      */
if (brackt && bound)
    mq = x + 0.66 * (y - x);
    if (x < y)
        if (mq < newt)
            newt = mq;
        else
            newt = mq;
        end
    end
end
%     /* Return the new trial value. */
t = newt;
info=0;
return;
end


