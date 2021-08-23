function OUT=States_Vclamp(t,O,I,Oinf,TauO,DAinf,TauDA,OSi,V)
Il=OSi(t);
OUT=[(Oinf([V,Il])-O)./TauO([V,Il]);
    (DAinf([V,Il])-I)./TauDA([V,Il])];
end