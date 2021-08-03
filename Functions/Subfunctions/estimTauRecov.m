function Tau = estimTauRecov(TauDArecov,tVr,powerb,options)
Tau = zeros(size(TauDArecov));
for iestimTau=1:size(tVr,2)
    X0 = TauDArecov(iestimTau)*powerb;
    if isfield(options,'tend') && isfield(options,'dt')
        to=0:options.dt:options.tend;
    else
        to = 0:X0/100:X0*20;
    end
    if isfield(options,'lb')
        LB = options.lb;
    else
        LB = 0;
    end
    if isfield(options,'ub')
        UB = options.ub;
    else
        UB = X0*20;
    end
    ty=(1-exp(-to./TauDArecov(iestimTau))).^powerb;
    OBJo=@(p2,t) (1-exp(-t/p2));
    opt=optimoptions(options.optimoptions{:});
    problemo=createOptimProblem('lsqcurvefit',...
        'objective',OBJo,...
        'xdata',to,'ydata',ty,...
        'x0',X0,...
        'lb',LB,...
        'ub',UB, 'options',opt);
    mso = MultiStart(options.msopt{:});
    [Tau(iestimTau)]= run(mso,problemo,options.msnr);
end
end