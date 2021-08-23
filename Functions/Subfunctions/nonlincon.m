function [c,ceq] = nonlincon(p,pa,pb,Oinffun,PidxO,DAinffun,PidxDA,Grectfun,PidxGr,XO,XO0,XDA,XDA0,V1,V2,extrac,extraceq,TauOoff,TauDAoff,ITau)

c=[  Oinffun(p(PidxO),XO)-1;
    -Oinffun(p(PidxO),XO);
    DAinffun(p(PidxDA),XDA)-1;
    -DAinffun(p(PidxDA),XDA);
    -Grectfun(p(PidxGr),V1);
    Grectfun(p(PidxGr),V2);
    1-DAinffun(p(PidxDA),ITau)-pa.*TauDAoff./(pa.*TauDAoff+pb.*TauOoff)]; %% this to make sure it doesnt go back up after tauoff

if extrac.IO_smaller_flag
    c = [c;Oinffun(p(PidxO),extrac.IO_smaller).^pa-extrac.Olim_smaller];
end
if extrac.IDA_smaller_flag
    c = [c;DAinffun(p(PidxDA),extrac.IDA_smaller).^pb-extrac.DAlim_smaller];
end
if extrac.IO_bigger_flag
    c = [c;-Oinffun(p(PidxO),extrac.IO_bigger).^pa+extrac.Olim_bigger];
end
if extrac.IDA_bigger_flag
    c = [c;-DAinffun(p(PidxDA),extrac.IDA_bigger).^pb+extrac.DAlim_bigger];
end

ceq=[DAinffun(p(PidxDA),XDA0)-1;
    Oinffun(p(PidxO),XO0)
    Grectfun(p(PidxGr),extraceq.VG)];

if any([isnan(c);isnan(ceq)])
    c(isnan(c)) = 1;
    ceq(isnan(ceq))=1;
end
if any([isinf(c);isinf(ceq)])
    c(isinf(c)) = 1;
    ceq(isinf(ceq))=1;
end
if any([imag(c)~=0;imag(ceq)~=0])
    c(imag(c)~=0) = 1;
    ceq(imag(ceq)~=0)=1;
end
end