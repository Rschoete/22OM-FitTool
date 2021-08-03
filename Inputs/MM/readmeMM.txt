modToff: Tau Off modified easier to satisfy nonlin con => tau off/4
modToff1711: change fit order GODAtDAtO
mdToff1914: also GODAtDAtO and different weighting fopt (not very succesful)
modToff2012: Also GODAtDAtO an [1000,1000,1000,1/25,1,250,10] weights + new scale system (fraction of init boundary)
MM_sga_t200416_modToff_d20042214: similar above but now ga
MM_spsoBC_t200416_modToff_d20042215&17: similar above but now psoBC difference probably nr generations


MM_spsoBC_t200416_modToff_20042921: max gen 100 new run with hybrid
MM_sga_t200416_modToff_20042921: max gen 100 new run with hybrid


above have
options_fopt.extrac_nonlcon.IO_smaller = [4000];
options_fopt.extrac_nonlcon.Olim_smaller = [0.4];
options_fopt.extrac_nonlcon.IO_bigger = [100];
options_fopt.extrac_nonlcon.Olim_bigger = [0.1];

below have
options_fopt.extrac_nonlcon.IO_smaller = [4000];
options_fopt.extrac_nonlcon.Olim_smaller = [0.5];
options_fopt.extrac_nonlcon.IO_bigger = [];
options_fopt.extrac_nonlcon.Olim_bigger = [];

MM_sga_t200416_modToff_20043015
MM_sms_t200416_modToff_20043015
MM_spsoBC_t200416_modToff_20043015


modToff2 all tau off eaual to 0.005
weights final opt 10000,1000,1000,1/25,1,250,10

d20051116 old weights  1000,1000,1000,1/25,1,250,10

20071714: same as above but now rerun with update fittool (startpoints in ga and psoBC and taurec updated)
Maxtime 8uur (run op wicasim!!! 20c vs 4 => compensate)
tOnL = > weights final opt 10000,1000,1000,1/25,1,250,10  (wTau on *10)

MM_sms_t200416_modToff2_d20072913_tOnL_GODAtDAtO  wTauon *10 order GODAtDAtO  GODA fixed ms tauda 2000 tauO2000 goda 3000 *uur max time run on wicasims!!
