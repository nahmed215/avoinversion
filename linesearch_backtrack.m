function [ x, g, f, info, count, img] = linesearch_backtrack(fg, xp,  f, g, s, stp, ftol,wolfe, max_linesearch, min_step, max_step,linesearch_type)

x=xp;
ret=0;
width=0;
dg=0;
norm=0;
finit=0;
dginit=0;
dgtest=0;
dec=.5;
inc=2.1;
count=0;
info=0;


%     /* Check the input parameters for errors. */
if (stp <= 0.)
    info=999;
    return;
end

%     /* Compute the initial gradient in the search direction. */
dginit=sum(g.*s);

%     /* Make sure that s points to a descent direction. */
if (0 < dginit)
    info=998;
    return;
end


%     /* The initial value of the objective function. */
finit = f;
dgtest = ftol * dginit;

while (1)
    x=xp;
    x=x+ stp*s;
    
    if (nargout>5)
        [f,g,img]=feval(fg,x);
    else
        [f,g]=feval(fg,x);
    end
    
    count=count+1;

    if (f > finit + stp * dgtest)
        width = dec;
    else
%         display('testing linesearch');
        %             /* The sufficient decrease condition (Armijo condition). */
        if (linesearch_type == 1)
            %                 /* Exit with the Armijo condition. */
%             display('Armijo');
            return;
        end
        
        
        % 	        /* Check the Wolfe condition. */
        dg=sum(g.*s);
        if (dg < wolfe * dginit)
            width = inc;
        else
            if(linesearch_type == 2)
                % 		            /* Exit with the regular Wolfe condition. */
%                 display('Wolfe');
                return;
            end
        end
        
        % 		        /* Check the strong Wolfe condition. */
        if(dg > -wolfe * dginit)
            width = dec;
        else
            % 		            /* Exit with the strong Wolfe condition. */
%             display('Strong wolfe');
            return;
        end
    end
    
    
    if (stp < min_step)
        %             /* The step is the minimum value. */
        info=997;
        return;
    end
    if (stp > max_step)
        %             /* The step is the maximum value. */
        info=996;
        return;
    end
    if (max_linesearch <= count)
        %             /* Maximum number of iteration. */
        info=995;
        return;
    end
    
    stp = stp*width;
end
end



