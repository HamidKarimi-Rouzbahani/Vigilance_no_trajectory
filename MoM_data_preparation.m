%% Sample response processing
clc
clear all;
subjects=1:20;
Blocks={'1_','2_','3_','4_','5_','6_','7_','8_','9_','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};
load('Order_of_blocks_for_Psych.mat');

Num_double_press=0;
Num_before_appearance_association=0;


for Subj_Num=[(subjects)]

    address='E:\Hamid Projects\Eye-tracking_lab\Vigilance\Newest_no_trajectory';
    dirs=dir(address);
    cd(address);

    blk_nums=30;


    for blk=1:blk_nums
        %         clearvars -except blk_nums Distances_traj Num_before_appearance_association Num_double_press sums nums BlockNumber Response_to_dots Cued_color_block Active_monitoring_block green_red_color up_down_direct_defl_det left_right_direct_app_det left_right_direct_app First_color_blocks Num_moving_dots Number_of_trials_in_blocks dprime_pattern_active dprime_pattern_monit TNR_active FNR_active FPR_active FPR_monit FPR_top FPR_bottom TPR_active TPR_monit TPR_top TPR_bottom correct_reaction_monit correct_reaction_active correct_reaction_times_top correct_reaction_times_bottom accuracy_active accuracy_monit accuracy_top accuracy_bottom Subj_Num dirs Blocks blk accuracy_top accuracy_bottom error_top error_bottom
        correct_reaction_times_att=0;
        for i=1:size(dirs,1)
            if Subj_Num<10
                num_chars=13;
            else
                num_chars=14;
            end
            if strncmp(dirs(i).name,['Subj_',num2str(Subj_Num),'_Blk_',Blocks{blk}],num_chars)
                load(dirs(i).name);
                break;
            end
        end
        ActMon=str2double(dirs(i).name(end-8));
        mean_sampling_time=1./60;
        for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
            tr=ceil(dot_num./Num_moving_dots);
            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;

            if ~isempty(find(key_pressed1(dot_in_trial,:,tr),1))

                if length(find(key_pressed1(dot_in_trial,:,tr)))>1
                    Num_double_press=Num_double_press+1;
                end

                key_press_sample=find(key_pressed1(dot_in_trial,:,tr), 1, 'first');
                if isnan(distance_traj1(dot_num,key_press_sample))
                    Num_before_appearance_association=Num_before_appearance_association+1;
                    distance_traj1(dot_num,key_press_sample)=3000;
                end
                dist_relative_to_boundary(dot_in_trial,tr)=distance_traj1(dot_num,key_press_sample)-hitting_border_distance;
            else
                dist_relative_to_boundary(dot_in_trial,tr)=nan;
            end
            distance_change_per_sample(dot_in_trial,tr)=(distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+10)-distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+20))./(11);

            if ~isempty(find(key_pressed2(dot_in_trial,:,tr),1))

                if length(find(key_pressed2(dot_in_trial,:,tr)))>1
                    Num_double_press=Num_double_press+1;
                end

                key_press_sample2=find(key_pressed2(dot_in_trial,:,tr), 1, 'first' );
                if isnan(distance_traj2(dot_num,key_press_sample2))
                    Num_before_appearance_association=Num_before_appearance_association+1;
                    distance_traj2(dot_num,key_press_sample)=3000;
                end
                dist_relative_to_boundary2(dot_in_trial,tr)=distance_traj2(dot_num,key_press_sample2)-hitting_border_distance;
            else
                dist_relative_to_boundary2(dot_in_trial,tr)=nan;
            end
            distance_change_per_sample2(dot_in_trial,tr)=(distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+10)-distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+20))./(11);
        end


        distance_change_per_sample(distance_change_per_sample<0)=mean(distance_change_per_sample(distance_change_per_sample>0));
        distance_change_per_sample2(distance_change_per_sample2<0)=mean(distance_change_per_sample2(distance_change_per_sample2>0));

        reaction_times=((-dist_relative_to_boundary)./distance_change_per_sample).*mean_sampling_time;
        reaction_times2=((-dist_relative_to_boundary2)./distance_change_per_sample2).*mean_sampling_time;
        %% Behavioural Performance

        % attended
        tp_att=0;
        tn_att=0;
        fp_F_att=0;
        fp_S_att=0;
        fp_T_att=0;
        fn_att=0;

        g=0;
        for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
            tr=ceil(dot_num./Num_moving_dots);
            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;

            if sum(dot_in_trial==top_events(:,tr))==1 && dot_color(dot_in_trial,tr)==First_color_blocks(Subj_Num,blk)
                g=g+1;
                if isnan(reaction_times(dot_in_trial,tr)) && (top_events(tr)~=top_targets(tr))
                    tn_att=tn_att+1;    % number of non-target events with no resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                    fp_F_att=fp_F_att+1;    % number of non-target events with fast resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)>=0)
                    fp_S_att=fp_S_att+1;    % number of non-target events with Slow resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                    fp_T_att=fp_T_att+1;    % number of target events with Too early resp;
                elseif isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr)
                    fn_att=fn_att+1;    % number of target events with no resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && reaction_times(dot_in_trial,tr)>0
                    tp_att=tp_att+1;    % number of target events with resp;
                    correct_reaction_times_att=correct_reaction_times_att+reaction_times(dot_in_trial,tr);
                end
            end

            if sum(dot_in_trial==top_events2(:,tr))==1 && dot_color2(dot_in_trial,tr)==First_color_blocks(Subj_Num,blk)
                g=g+1;
                if isnan(reaction_times2(dot_in_trial,tr)) && (top_events2(tr)~=top_targets2(tr))
                    tn_att=tn_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)
                    fp_F_att=fp_F_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)>=0)
                    fp_S_att=fp_S_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)% || reaction_times2(dot_in_trial,tr)>time_to_touch_the_obstacle2(dot_in_trial,tr))
                    fp_T_att=fp_T_att+1;
                elseif isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr)
                    fn_att=fn_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && reaction_times2(dot_in_trial,tr)>0 %&& reaction_times2(dot_in_trial,tr)<time_to_touch_the_obstacle2(dot_in_trial,tr)
                    tp_att=tp_att+1;
                    correct_reaction_times_att=correct_reaction_times_att+reaction_times2(dot_in_trial,tr);
                end
            end
        end

        fp_att=fp_F_att+fp_S_att+fp_T_att;

        correct_reaction_times_att=correct_reaction_times_att./tp_att;



        % Removed the unattended dots for simplicity of the data
        blk_all=blk;
        if ActMon==1
            cond=1;
        else
            cond=2;
        end
        Data{cond,1}(blk_all,Subj_Num)=(tp_att)./(tp_att+fn_att);
        % False alarm
        Data{cond,2}(blk_all,Subj_Num)=(fp_att)./(fp_att+tn_att);
        % Reaction time
        Data{cond,3}(blk_all,Subj_Num)=correct_reaction_times_att;
        if ActMon==1
            cond=2;
            Data{cond,1}(blk_all,Subj_Num)=nan;
            Data{cond,2}(blk_all,Subj_Num)=nan;
            Data{cond,3}(blk_all,Subj_Num)=nan;
        else
            cond=1;
            Data{cond,1}(blk_all,Subj_Num)=nan;
            Data{cond,2}(blk_all,Subj_Num)=nan;
            Data{cond,3}(blk_all,Subj_Num)=nan;            
        end
    end
    [Subj_Num]
end


%% Saving data as Excel file for analysis
for Subj=subjects
    for cond=1:size(Data,1)
        if ~ismember(Subj,subjects)
            Data{cond,1}(:,Subj)=nan;
            Data{cond,2}(:,Subj)=nan;
            Data{cond,3}(:,Subj)=nan;
        end
    end
    cond=1;
    Hit_rate_condition=Data{cond,1}(:,Subj); % Hit rate in condition
    Mean_Hit_rate=nanmean(Hit_rate_condition);
    FA_rate_condition=Data{cond,2}(:,Subj); % FA rate in condition
    Mean_FA_rate=nanmean(FA_rate_condition);
    Reaction_time_condition=Data{cond,3}(:,Subj); % reaction time in condition
    Mean_Reaction_time=nanmean(Reaction_time_condition);
    
    HR=[Hit_rate_condition;nan(5,1);Mean_Hit_rate];
    FA=[FA_rate_condition;nan(5,1);Mean_FA_rate];
    RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];
    T = table(HR,FA,RT);
    T.Properties.VariableNames = {['Hit_rate_Active'] ['FA_rate_Active'] ['Reaction_Time_Active']};
    Ttotal=T;
    Data_csv_total=[HR FA RT];
    
    for cond=2:size(Data,1)
        
        Hit_rate_condition=Data{cond,1}(:,Subj); % Hit rate in condition
        Mean_Hit_rate=nanmean(Hit_rate_condition);
        FA_rate_condition=Data{cond,2}(:,Subj); % FA rate in condition
        Mean_FA_rate=nanmean(FA_rate_condition);
        Reaction_time_condition=Data{cond,3}(:,Subj); % reaction time in condition
        Mean_Reaction_time=nanmean(Reaction_time_condition);

        HR=[Hit_rate_condition;nan(5,1);Mean_Hit_rate];
        FA=[FA_rate_condition;nan(5,1);Mean_FA_rate];
        RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];
        T = table(HR,FA,RT);
        T.Properties.VariableNames = {['Hit_rate_Monit'] ['FA_rate_Monit'] ['Reaction_Time_Monit']};

        Ttotal=[Ttotal T];
        Data_csv_total=horzcat(Data_csv_total,[HR RT]);
    end

    filename = ['MoM_data_no_trajectory.xlsx']; % Change the name to anything you prefer
    writetable(Ttotal,filename,'Sheet',['Subj_' num2str(Subj)])
       
    [Subj]
end

%% Plotting the results

plott=1; % 1= Hit Rate; 2= False Alarm; 3=RT
gca = axes('Position',[0.22 0.25 0.775 0.72]);
if plott==1
    data_to_plot1=(Data{1,plott}(:,subjects))*100;
    data_to_plot2=(Data{2,plott}(:,subjects))*100;
    miny=20;
    maxy=100;
    yticks=[0:20:100];
elseif plott==2
    data_to_plot1=(Data{1,plott}(:,subjects))*100;
    data_to_plot2=(Data{2,plott}(:,subjects))*100;
    miny=0;
    maxy=40;
    yticks=[0:20:40];
elseif plott==3
    data_to_plot1=Data{1,3}(:,subjects)*1000;
    data_to_plot2=Data{2,3}(:,subjects)*1000;
    miny=200;
    maxy=400;
    yticks=[200:50:400];
end
MeanA=(nanmean([data_to_plot1(1:15,:) data_to_plot1(16:30,:)],2));
StdA=nanstd([data_to_plot1(1:15,:) data_to_plot1(16:30,:)]');
MeanM=(nanmean([data_to_plot2(1:15,:) data_to_plot2(16:30,:)],2));
StdM=nanstd([data_to_plot2(1:15,:) data_to_plot2(16:30,:)]');
plots(1)=plot([1:15],MeanA');
hold on;
plots(2)=plot([1:15],MeanM');

xticks={'','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',''};
xlim([0 16])
ylim([miny maxy])
if plott==1
    ylabel({'Hit rate (%)'})
elseif plott==2
    ylabel({'False alarm (%)'})
elseif plott==3
    ylabel({'Reaction time (ms)'})
end
set(gca,'FontSize',20,'LineWidth',3,'XTick',...
    [0:16],'XTickLabel',{'','','','','','','','','','','','','','','','',''},...
    'YTick',yticks,'YTickLabel',num2str(yticks'),'ycolor','k','xcolor','k');
box off;
if plott==1
    legend ({'Active','Monitoring'},'location','southwest','edgecolor','none')
end
