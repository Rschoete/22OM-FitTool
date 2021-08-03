function [state,options,optchanged] = stopPSO(options,state,flag)
stop=false;
msg=lastwarn;
optchanged = 0;
if contains(msg,'singular') || contains(msg,'Unable to meet integration tolerances') || contains(msg,'adversely affect')
    stop=true;
    state.StopFlag='y';
    options.StallGenLimit=0;
    options.StallTimeLimit=0;
    options.TimeLimit=0;
    lastwarn('');
    optchanged = 1;
end
end
