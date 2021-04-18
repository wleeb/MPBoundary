        function main
        prini(13,1);
%
        randn(1,1000);

        m=1000;
        n=2000;

        prinf('m=',m,1);
        prinf('n=',n,1);

        gam = m / n;

        as = rand(1,m);
        bs = rand(1,n);

        awhts=ones(1,m);
        awhts = awhts / sum(awhts);
        bwhts=ones(1,n);
        bwhts = bwhts / sum(bwhts);


        difa = sum(awhts) - 1;
        difb = sum(bwhts) - 1;
        prin2('difa=',difa,1);
        prin2('difb=',difb,1);

%%%        prinstop;

   
        bedge = mpbdry_edge(as,bs,awhts,bwhts,m,n,gam);

        prin2('bedge=',bedge,1);


        krank=3;


        ells = bedge + rand(1,krank);
        ells = sort(ells,'descend');
        [y,x,ep,u,v] = genspiket(ells,m,n,krank,as,awhts,bs,bwhts);


        [uy,sy,vy] = svshr_svdsmart(y,m,n,krank);


        chk0 = norm(y*vy - uy*diag(sy),'fro') / norm(sy,'fro');

        prin2('chk0=',chk0,1);



%
%        compare asyptotic and observed edge
%
        bedge_obs = norm(ep)^2;
        err = (bedge_obs - bedge) / bedge;

        prin2('edge,observed=',bedge_obs,1);
        prin2('edge,exact=',bedge,1);
        prin2('relative error=',err,1);




%
%        compare spiked model parameters
%
        i=2;
        [rlam,cout,cinn] = mpbdry_sforw(ells(i),as,bs,...
            awhts,bwhts,m,n,gam);
        rlam_obs = sy(i)^2;
        err = (rlam_obs - rlam) / rlam;

        prin2('lambda,observed=',rlam_obs,1);
        prin2('lambda,exact=',rlam,1);
        prin2('relative error=',err,1);


%%%        prinstop;


        cout_obs = abs(sum(u(:,i).*uy(:,i)));
        err = (cout_obs - cout) / cout;


        prin2('cosine,observed=',cout_obs,1);
        prin2('cosine,exact=',cout,1);
        prin2('relative error=',err,1);


        cinn_obs = abs(sum(v(:,i).*vy(:,i)));
        err = (cout_obs - cout) / cout;


        prin2('cosine,observed=',cinn_obs,1);
        prin2('cosine,exact=',cinn,1);
        prin2('relative error=',err,1);

%%%        prinstop;

%
%        check spike backwards mapping
%
        [ell2,cout2,cinn2] = mpbdry_sback(rlam,as,bs,...
            awhts,bwhts,m,n,gam);

        chk0 = abs(ell2 - ells(i));
        chk0 = chk0+abs(cout2 - cout);
        chk0 = chk0+abs(cinn2 - cinn);


        prin2('chk0, backwards mapping=',chk0,1);




        prinstop;
        end
%
%
%
%
%
        function [y,x,ep,u,v] = genspiket(ells,m,n,k,as,awhts,bs,bwhts)
%
        ep = randn(m,n) / sqrt(n);
        ep = diag(sqrt(as)) * ep * diag(sqrt(bs));

        u = randn(m,k);
        v = randn(n,k);


        u = gramschmidt(u,m,k);
        v = gramschmidt(v,n,k);

        chk0 = norm(u'*u - eye(k),'fro');
        chk0 = chk0 + norm(v'*v - eye(k),'fro');

        prin2('chk0=',chk0,1);
%%%        prinstop;


        x = u * diag(sqrt(ells)) * v';

        xyz = svds(x,k)';
        prin2('xyz=',xyz,k);
        prin2('xyz=',sqrt(ells),k);
%%%        prinstop;

        y = x + ep; 

        end
%
%
%
%
%
        function u = gramschmidt(uin,m,k)
%
        u = uin;

        u(:,1) = u(:,1) / norm(u(:,1));

%%%        dif = norm(u(:,1))-1;
%%%        prin2('dif=',dif,1);
%%%        prinstop;

        for i=2:k
%
        for ijk = 1:2
%
        for j=1:i-1
%
        pij = sum(u(:,i).*u(:,j));
        u(:,i) = u(:,i) - pij*u(:,j);
    end
        u(:,i) = u(:,i) / norm(u(:,i));
    end

    end

        chk0 = norm(u'*u - eye(k),'fro');
        prin2('chk0=',chk0,1);
%%%        prinstop;

        end
%
%
%
%
%
        function stra = sbar2stra(sbar,z,gam)
%
%        given sbar(z), computes stra(z) (Stieltjes transform)
%
        stra = (sbar + (1-gam)/z) / gam;

        end
%
%
%
%
%
        function stra = stieltjes_mp(z,gam,sig)
%
%        evaluate the Stieltjes transform of standard MP law
%
        b = sig^2*(1-gam) - z;
        disc = (z - sig^2 - gam*sig^2)^2 - 4*gam*sig^4;
        a = 2*gam*z*sig^2;

        stra = (b + sqrt(disc)) / a;

        end
%
%
%
%
%
        function [stra,stra_der] = svals2integrs(s,m,n,rlam);
%
%        computes empirical estimates of the stieljes and D-transforms 
%        at the value rlam, with empirical singular values s (already
%        normalized by dimension).
%
        stra = 0;
        stra_der = 0;

        lmin = min(m,n);
%
        for i=1:lmin
%
        stra = stra + 1 / (s(i)^2 - rlam);
        stra_der = stra_der + 1 / (rlam - s(i)^2)^2;        
    end

        stra = stra/lmin;
        stra_der = stra_der/lmin;

        if (n <= m)
%
        stra = (n/m)*stra - (m-n)/m/rlam;
        stra_der = (n/m)*stra_der + (m-n)/m/rlam^2;
    end

        end
%
%
%
%
%
        function [u,s,v] = svshr_svdsmart(a,m,n,k)
%
        if (m/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
    end

        if (m/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
    end

        end
%
%
%
%
%
