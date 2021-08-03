function out = adjFields(original,new,level)
%function that merges structures
fn_lvl1_org = fieldnames(original);
fn_lvl1_new = fieldnames(new);
idx = 0;
newfields = {};


out = new;
for ifnn = 1:length(fn_lvl1_new)
    if ~any(strcmp(fn_lvl1_org,fn_lvl1_new{ifnn}))
        idx = idx+1;
        newfields{idx} =  fn_lvl1_new{ifnn};
    end
end

if idx>0
    fprintf('\n new fields added at lvl 1: %s \n',strjoin(newfields,', '))
    newfields = {};
    idx=0;
end

for ifno = 1:length(fn_lvl1_org)
    if ~any(strcmp(fn_lvl1_new,fn_lvl1_org{ifno}))
        out.(fn_lvl1_org{ifno}) = original.(fn_lvl1_org{ifno});
    else
        if level>1
            if isstruct(new.(fn_lvl1_org{ifno}))
                if ~isstruct(original.(fn_lvl1_org{ifno}))
                    error('mismatch orignal and new lvl1')
                end
                fn_lvl2_new = fieldnames(new.(fn_lvl1_org{ifno}));
                fn_lvl2_org = fieldnames(original.(fn_lvl1_org{ifno}));

                for ifnn = 1:length(fn_lvl2_new)
                    if ~any(strcmp(fn_lvl2_org,fn_lvl2_new{ifnn}))
                        idx = idx+1;
                        newfields{idx} =  fn_lvl2_new{ifnn};
                    end
                end

                if idx>0
                    fprintf('\n new fields added at lv2 1: %s \n',strjoin(newfields,', '))
                    newfields = {};
                    idx = 0;
                end

                for ifno_lvl2 = 1:length(fn_lvl2_org)
                    if ~any(strcmp(fn_lvl2_new,fn_lvl2_org{ifno_lvl2}))
                        out.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}) = original.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2});
                    else
                        if level>2
                            if isstruct(new.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}))
                                if ~isstruct(original.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}))
                                    error('mismatch original new lvl2')
                                end
                                fn_lvl3_new = fieldnames(new.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}));
                                fn_lvl3_org = fieldnames(original.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}));

                                for ifnn = 1:length(fn_lvl3_new)
                                    if ~any(strcmp(fn_lvl3_org,fn_lvl3_new{ifnn}))
                                        idx = idx+1;
                                        newfields{idx} =  fn_lvl3_new{ifnn};
                                    end
                                end

                                if idx>0
                                    fprintf('\n new fields added at lvl3 : %s \n',strjoin(newfields,', '))
                                    newfields = {};
                                    idx = 0;
                                end

                                for ifno_lvl3 = 1:length(fn_lvl3_org)
                                    if ~any(strcmp(fn_lvl3_new,fn_lvl3_org{ifno_lvl3}))
                                        out.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}).(fn_lvl3_org{ifno_lvl3}) = original.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}).(fn_lvl3_org{ifno_lvl3});
                                    else
                                        if level>3
                                            error('not encoded')
                                        end
                                    end

                                end
                            else
                                if isstruct(original.(fn_lvl1_org{ifno}).(fn_lvl2_org{ifno_lvl2}))
                                    error('mismatch original new lvl2')
                                end
                            end
                        end
                    end

                end
            else
                if isstruct(original.(fn_lvl1_org{ifno}))
                    error('mismatch orignal and new lvl1')
                end
            end
        end
    end
end