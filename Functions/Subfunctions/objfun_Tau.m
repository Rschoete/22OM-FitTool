function Out = objfun_Tau(p,xdata_Ilnz,ydata_Ilnz,xdata_Ilz,ydata_Ilz,wIlnz,wIlz,Taufun)
%Ilnz light intesity not zero
%Ilz light intesity zero
%used in fitTau
OUT = [wIlnz.*((ydata_Ilnz-Taufun(p,xdata_Ilnz))./ydata_Ilnz).^2;
    wIlz.*((ydata_Ilz-Taufun(p,xdata_Ilz))./ydata_Ilz).^2];
Out = sum(OUT);
end
