        function [val,der] = mpbdry_fminder(rlam,as,bs,awhts,bwhts,m,n,gam)
%
%        Evaluates the minimum of value of F(lambda,e) (for fixed rlam),
%        and the derivative of the minimum as a function of lambda
%
        nsteps = 1000;
        evars = zeros(1,nsteps);
        fders = zeros(1,nsteps);
        fders2 = zeros(1,nsteps);

%
%        initialize to the left of the minimizer (where fder < 0)
%
        emin = mpbdry_leftend(rlam,as,bs,awhts,bwhts,m,n,gam);
        einit=emin+abs(emin) + 10;

        ifder2=0;
        for i=1:1000
%
        [fder,fder2] = mpbdry_evalfder2(einit,rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder2);

        if (fder > 0)
%
        einit = (einit + emin)/2;
    end
        if (fder <= 0)
%
        break;
    end
    end

%
%        find the minimum using Newton
%
        evars(1) = einit;
        ifder2=1;
        kstop=0;
        deps = mpbdry_machzero();
        tol=sqrt(deps)/100;
%
        for i=1:nsteps
%
        [fders(i),fders2(i)] = mpbdry_evalfder2(evars(i),rlam,as,bs,...
            awhts,bwhts,m,n,gam,ifder2);
        evars(i+1) = evars(i) - fders(i) / fders2(i);
%
        if (abs(fders(i)) < tol)
%
        kstop = kstop+1;
    end
        if (kstop == 2)
%
        nsteps = i;
        break;
    end
    end

        fders=fders(1:nsteps);
        fders2=fders2(1:nsteps);
        evars=evars(1:nsteps);

        val = evars(nsteps);
%
%        evaluate the derivative of the minimum, wrt lambda
%
        [fder2,fdmix] = mpbdry_evalfmix(val,rlam,as,bs,awhts,bwhts,m,n,gam);
        der = -fdmix / fder2;


        end
%
