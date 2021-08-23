function stop = outfun(x,optimValues,state,plot_flag,OinfIfun,pOIidx,VarOinfI,DAinfIfun,pDAIidx,VarDAinfI,TauO,TauDA,Grectfun,pGidx,VarGrect,Vm_valsp,Target_sp,options,powera,powerb,Colors,plotTauRecov_flag,Vm_valtr,Target_tr)
stop=false;
fignr = get(gcf,'number');
if plot_flag && options.plot_flag
    fun_intm.Oinf = @(X) OinfIfun(x(pOIidx),X(:,VarOinfI));
    fun_intm.DAinf = @(X) DAinfIfun(x(pDAIidx),X(:,VarDAinfI));
    fun_intm.TauO = @(X) TauO(x,X);
    fun_intm.TauDA = @(X) TauDA(x,X);
    fun_intm.Grect = @(X) Grectfun(x(pGidx),X(:,VarGrect));
    runODE(Vm_valsp,Target_sp,options.OdeOpts,fun_intm,powera,powerb,Colors,plotTauRecov_flag,Vm_valtr,Target_tr,'fignr',fignr)

    mtit(['iter: ',num2str(optimValues.iteration)])

end
msg=lastwarn;
if contains(msg,'singular') || contains(msg,'Unable to meet integration tolerances') || contains(msg,'adversely affect')
    stop=true;
    lastwarn('');
end

end