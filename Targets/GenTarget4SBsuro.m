% generate surrogate data ChR2H134R 4SB
close all
clear all
clc
addpath(pwd)
addpath(genpath('../Functions/Subfunctions'))
oldd = cd('../Functions');
%BUG function does not exist anymore located elsewhere
funhandle = 'ChR2H134R';
Vm_val = -80:20:40;
Il_val = logspace(1,4,7);
OSpstart = 0.5;
OSpdsp = 0.5;
intervals = [0.1,0.5,1,2,3,5,7,10,13,15];
tend_sp = 1.5;
tend_2p = 2*(intervals+OSpdsp)+OSpstart;
Input.Model = '4SB';
plot_flag = 1;

Target=gentTarget_surrogateData(funhandle,Vm_val,Il_val,OSpstart,OSpdsp,intervals,tend_sp,tend_2p,Input,plot_flag);

cd(oldd);

if input('save (1/0)')
    save(['Target4SBsuro_',datestr(now,'yymmdd'),'.mat'],'Target')
end