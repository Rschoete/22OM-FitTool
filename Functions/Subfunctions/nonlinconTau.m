function [c,ceq] = nonlinconTau(p,pa,pb,TauO,XTauO,DAinf,TauDAoff)

c=[1-DAinf-pa.*TauDAoff./(pb.*TauO(p,XTauO)+pa.*TauDAoff)]; %% this to make sure it doesnt go back up after tauoff 

ceq=[];

if any([isnan(c);isnan(ceq)])
    c(isnan(c)) = 1;
    ceq(isnan(ceq))=1;
end

end