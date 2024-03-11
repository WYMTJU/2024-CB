% This code is used to get the spike after position stretching
clear all
clc
close all
% limit the posang_ontrack with running speed > vel_threshold
vel_candidate = [5,0]; % cm/s
file_video = 'Data_video_turn.mat';
for mult =[1.5,2,2.5,3]
    for vel_threshold=vel_candidate
        for nbrgion = 1
            if nbrgion ==1
                TTList0= 'TTList_dCA1_pyr.txt';
                file_spike = 'Data_spikes_dCA1_turn.mat';
                region = 'dCA1';
            else
                TTList0= 'TTList_mPFC_pyr.txt';
                file_spike = 'Data_spikes_mPFC_turn.mat';
                region = 'mPFC';
            end
            for icon=1
                if icon==1
                    directories_allData_control
                    mark='Con';
                else
                    directories_allData_Alzheimer % AD;
                    mark='AD';
                end
                file_com = strcat('Cells_doall_ccw_ratemap_vel_',num2str(vel_threshold),'_',region ,'_turn.mat');
                file_out = ['Data_spikes_stretch_x_t_',num2str(mult),'.mat'];
                for ns = 1:isession
                    path_ns = path{ns};
                    cd(path_ns);
                    disp(path_ns)
                    trackdata_ns = trackdata{ns};
                    if nbrgion == 1
                        if isempty(CSClist_CA1{ns})
                            fprintf(['=====day %.3g do not have hippocampal cells=====' '\n'], ns)
                            continue
                        end
                    else
                        if isempty(CSClist_PFC{ns})
                            fprintf(['=====day %.3g do not have prefrontal cells=====' '\n'], ns)
                            continue
                        end
                    end
                    load(file_spike,'spikes');
                    load(file_video);
                    load(file_com,'fieldProp_allccw');
                    
                    load(trackdata_ns,'Ts_prerunning','Ts_test');
                    Ts_start_stop = Ts_prerunning{1,2}./ 1000000;
                    Ts_start_stop = [Ts_start_stop; Ts_test{1,1}(:,[1,3])./ 1000000];
                    video_t = data_video(:,1);
                    video_x = data_video(:,4);
                    spikes_stretch_x = spikes; 
                    spikes_stretch_x_t = spikes; 
                    
                    for nc = 1:size(spikes,1)
                        t = spikes{nc,2};
                        x = spikes{nc,3}(:,3);
                        if length(cell2mat(struct2cell(fieldProp_allccw{nc}(1))))<2 % 没有spike，没有特征，则跳过
                            continue
                        end
                        if  fieldProp_allccw{nc}(1).size>pi   %If the size is relatively large, do not stretch
                            spikes_stretch_x_t{nc,2} = t;
                            spikes_stretch_x_t{nc,3}(:,3) = x;
                            continue
                        end
                        for nl = 1:length(Ts_start_stop)
                            com = fieldProp_allccw{nc}(1).x_COM;
                            ind0 = find(t>=Ts_start_stop(nl,1) & t<=Ts_start_stop(nl,2));
                            if isempty(ind0)
                                continue
                            end
                            x_lap = x(ind0); t_lap = t(ind0);
                            diff = wrapToPi(x_lap-com);
                            ind1 = find(diff>0);
                            ind2 = find(diff<0);
                            x2=[];
                            x2(ind1) = x_lap(ind1) + (diff(ind1)) *mult;
                            x2(ind2) = x_lap(ind2) + (diff(ind2)) *mult;
                            x2 = wrapTo2Pi(x2); 
                            spikes_stretch_x{nc,3}(ind0,3) = x2;
                            spikes_stretch_x_t{nc,3}(ind0,3) = x2;
                            
                            video_ind = find(video_t>=Ts_start_stop(nl,1) & video_t<=Ts_start_stop(nl,2));
                            video_x_lap = video_x(video_ind);
                            video_t_lap = video_t(video_ind);
                            t2=[];
                            for i=1:length(x2)
                                [~,ind_min]=min(abs(video_x_lap-x2(i)));
                                t2(i,1) = video_t_lap(ind_min);
                            end
                            spikes_stretch_x_t{nc,2}(ind0) = t2;
                        end
                        close all
                    end
                    save(file_out,'spikes_stretch_x_t','spikes_stretch_x');
                end
            end
        end
    end
end
