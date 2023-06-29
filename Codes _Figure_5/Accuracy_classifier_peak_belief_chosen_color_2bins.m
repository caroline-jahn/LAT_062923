clear all;

fsroot='/Volumes/buschman';

event='target';

for i=1:27
    initial_window=-500;
    event_list{i}='target';
    window_start_list(i)=initial_window+(i-1)*50;
end

task='Learning_Attentional_Templates';
subtask='exploreexploit/Restricted_FEF/Reset_RW_model';

window_size=200;
subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

N_boot=100;

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

color_for_ROI=[0 0.447 0.941; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')

color_belief(1,:)=colors(1,:);
color_belief(2,:)=colors(14,:);
color_belief(3,:)=colors(27,:);

%LIP
ROI='LIP';

for n_tt=1:N_boot
    
    for this_time=1:length(window_start_list)
        
        results_name=sprintf('Pseudo_pop_chosen_color_peak_template_2bins_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);
        if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
            load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_chosen_color_each_belief_correct','Classification_chosen_color_correct')
            
            LIP.Classification_belief_correct(this_time,:,n_tt)=Classification_belief_correct(:);
            
            LIP.Classification_correct_chosen_color(this_time,:,n_tt)=Classification_chosen_color_correct(:);
            
            LIP.Classification_correct_chosen_color_each_belief(this_time,:,:,:,n_tt)=Classification_chosen_color_each_belief_correct(:,:,:);
            
            clear Classification*
        else
            sprintf('%s is missing',results_name)
            
            LIP.Classification_belief_correct(this_time,1:54,n_tt)=NaN;
            
            LIP.Classification_correct_chosen_color(this_time,1:54,n_tt)=NaN;
            
            LIP.Classification_correct_chosen_color_each_belief(this_time,1:18,1:3,1:3,n_tt)=NaN;
            
            
        end
    end
end

%FEF
ROI='FEF';

for n_tt=1:N_boot
    
    for this_time=1:length(window_start_list)
        
        results_name=sprintf('Pseudo_pop_chosen_color_peak_template_2bins_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);
        if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
            load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_chosen_color_each_belief_correct','Classification_chosen_color_correct')
            
            FEF.Classification_belief_correct(this_time,:,n_tt)=Classification_belief_correct(:);
            
            FEF.Classification_correct_chosen_color(this_time,:,n_tt)=Classification_chosen_color_correct(:);
            
            FEF.Classification_correct_chosen_color_each_belief(this_time,:,:,:,n_tt)=Classification_chosen_color_each_belief_correct(:,:,:);
            
            clear Classification*
        else
            sprintf('%s is missing',results_name)
            
            FEF.Classification_belief_correct(this_time,1:54,n_tt)=NaN;
            
            FEF.Classification_correct_chosen_color(this_time,1:54,n_tt)=NaN;
            
            FEF.Classification_correct_chosen_color_each_belief(this_time,1:18,1:3,1:3,n_tt)=NaN;
            
            
        end
    end
end

%PFC
ROI='PFC';

for n_tt=1:N_boot
    
    for this_time=1:length(window_start_list)
        
        results_name=sprintf('Pseudo_pop_chosen_color_peak_template_2bins_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);
        if exist(sprintf('%s.mat',fullfile(data_path_clasifier,results_name)))==2
            load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_chosen_color_each_belief_correct','Classification_chosen_color_correct')
            
            PFC.Classification_belief_correct(this_time,:,n_tt)=Classification_belief_correct(:);
            
            PFC.Classification_correct_chosen_color(this_time,:,n_tt)=Classification_chosen_color_correct(:);
            
            PFC.Classification_correct_chosen_color_each_belief(this_time,:,:,:,n_tt)=Classification_chosen_color_each_belief_correct(:,:,:);
            
            clear Classification*
        else
            sprintf('%s is missing',results_name)
            
            PFC.Classification_belief_correct(this_time,1:54,n_tt)=NaN;
            
            PFC.Classification_correct_chosen_color(this_time,1:54,n_tt)=NaN;
            
            PFC.Classification_correct_chosen_color_each_belief(this_time,1:18,1:3,1:3,n_tt)=NaN;
            
            
        end
    end
end


%%

for i=1:length(window_start_list)
    p_belief_LIP(i)=z_test_function_bootstrap(nanmean(LIP.Classification_belief_correct(i,:,:),2),1/2);
    p_belief_FEF(i)=z_test_function_bootstrap(nanmean(FEF.Classification_belief_correct(i,:,:),2),1/2);
    p_belief_PFC(i)=z_test_function_bootstrap(nanmean(PFC.Classification_belief_correct(i,:,:),2),1/2);
    
    p_chosen_color_LIP(i)=z_test_function_bootstrap(nanmean(LIP.Classification_correct_chosen_color(i,:,:),2),1/2);
    p_chosen_color_FEF(i)=z_test_function_bootstrap(nanmean(FEF.Classification_correct_chosen_color(i,:,:),2),1/2);
    p_chosen_color_PFC(i)=z_test_function_bootstrap(nanmean(PFC.Classification_correct_chosen_color(i,:,:),2),1/2);
    for l=1:2
        for m=1:2
            p_chosen_color_each_belief_LIP(i,l,m)=z_test_function_bootstrap(nanmean(LIP.Classification_correct_chosen_color_each_belief(i,:,l,m,:),2),1/2);
            p_chosen_color_each_belief_FEF(i,l,m)=z_test_function_bootstrap(nanmean(FEF.Classification_correct_chosen_color_each_belief(i,:,l,m,:),2),1/2);
            p_chosen_color_each_belief_PFC(i,l,m)=z_test_function_bootstrap(nanmean(PFC.Classification_correct_chosen_color_each_belief(i,:,l,m,:),2),1/2);
        end
    end
end


%%
for i=1:3
    color_for_ROI_early(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_early(i,:)=hsv2rgb(color_for_ROI_early(i,1),color_for_ROI_early(i,2),1);
    color_for_ROI_late(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_late(i,:)=hsv2rgb(color_for_ROI_late(i,1),color_for_ROI_late(i,2),0.55);
    
end

%%
for n_tt=1:N_boot
    for this_time=1:length(window_start_list)
        for m=1:2
            LIP.Classification_correct_chosen_color_same_belief(this_time,1+(m-1)*size(PFC.Classification_correct_chosen_color_each_belief,2):m*size(PFC.Classification_correct_chosen_color_each_belief,2),n_tt)=LIP.Classification_correct_chosen_color_each_belief(this_time,:,m,m,n_tt);
            FEF.Classification_correct_chosen_color_same_belief(this_time,1+(m-1)*size(PFC.Classification_correct_chosen_color_each_belief,2):m*size(PFC.Classification_correct_chosen_color_each_belief,2),n_tt)=FEF.Classification_correct_chosen_color_each_belief(this_time,:,m,m,n_tt);
            PFC.Classification_correct_chosen_color_same_belief(this_time,1+(m-1)*size(PFC.Classification_correct_chosen_color_each_belief,2):m*size(PFC.Classification_correct_chosen_color_each_belief,2),n_tt)=PFC.Classification_correct_chosen_color_each_belief(this_time,:,m,m,n_tt);
        end
    end
end

%%
count=1;

for m=1:2
    for n=1:2
        if m~=n
            for this_time=1:length(window_start_list)
                for n_tt=1:N_boot
                    LIP.Classification_correct_chosen_color_diff_belief(this_time,1+(count-1)*size(PFC.Classification_correct_chosen_color_each_belief,2):count*size(PFC.Classification_correct_chosen_color_each_belief,2),n_tt)=LIP.Classification_correct_chosen_color_each_belief(this_time,:,m,n,n_tt);
                    FEF.Classification_correct_chosen_color_diff_belief(this_time,1+(count-1)*size(PFC.Classification_correct_chosen_color_each_belief,2):count*size(PFC.Classification_correct_chosen_color_each_belief,2),n_tt)=FEF.Classification_correct_chosen_color_each_belief(this_time,:,m,n,n_tt);
                    PFC.Classification_correct_chosen_color_diff_belief(this_time,1+(count-1)*size(PFC.Classification_correct_chosen_color_each_belief,2):count*size(PFC.Classification_correct_chosen_color_each_belief,2),n_tt)=PFC.Classification_correct_chosen_color_each_belief(this_time,:,m,n,n_tt);
                end
            end
            count=count+1;
        end
    end
end



%%


LIP.Same_belief_chosen_color=horzcat(reshape(LIP.Classification_correct_chosen_color_each_belief(:,:,1,1,:),size(LIP.Classification_correct_chosen_color_each_belief,1),size(LIP.Classification_correct_chosen_color_each_belief,2),size(LIP.Classification_correct_chosen_color_each_belief,5)),...
    reshape(LIP.Classification_correct_chosen_color_each_belief(:,:,2,2,:),size(LIP.Classification_correct_chosen_color_each_belief,1),size(LIP.Classification_correct_chosen_color_each_belief,2),size(LIP.Classification_correct_chosen_color_each_belief,5)));

LIP.Other_belief_chosen_color=horzcat(reshape(LIP.Classification_correct_chosen_color_each_belief(:,:,1,2,:),size(LIP.Classification_correct_chosen_color_each_belief,1),size(LIP.Classification_correct_chosen_color_each_belief,2),size(LIP.Classification_correct_chosen_color_each_belief,5)),...
    reshape(LIP.Classification_correct_chosen_color_each_belief(:,:,2,1,:),size(LIP.Classification_correct_chosen_color_each_belief,1),size(LIP.Classification_correct_chosen_color_each_belief,2),size(LIP.Classification_correct_chosen_color_each_belief,5)));

FEF.Same_belief_chosen_color=horzcat(reshape(FEF.Classification_correct_chosen_color_each_belief(:,:,1,1,:),size(FEF.Classification_correct_chosen_color_each_belief,1),size(FEF.Classification_correct_chosen_color_each_belief,2),size(FEF.Classification_correct_chosen_color_each_belief,5)),...
    reshape(FEF.Classification_correct_chosen_color_each_belief(:,:,2,2,:),size(FEF.Classification_correct_chosen_color_each_belief,1),size(FEF.Classification_correct_chosen_color_each_belief,2),size(FEF.Classification_correct_chosen_color_each_belief,5)));

FEF.Other_belief_chosen_color=horzcat(reshape(FEF.Classification_correct_chosen_color_each_belief(:,:,1,2,:),size(FEF.Classification_correct_chosen_color_each_belief,1),size(FEF.Classification_correct_chosen_color_each_belief,2),size(FEF.Classification_correct_chosen_color_each_belief,5)),...
    reshape(FEF.Classification_correct_chosen_color_each_belief(:,:,2,1,:),size(FEF.Classification_correct_chosen_color_each_belief,1),size(FEF.Classification_correct_chosen_color_each_belief,2),size(FEF.Classification_correct_chosen_color_each_belief,5)));

PFC.Same_belief_chosen_color=horzcat(reshape(PFC.Classification_correct_chosen_color_each_belief(:,:,1,1,:),size(PFC.Classification_correct_chosen_color_each_belief,1),size(PFC.Classification_correct_chosen_color_each_belief,2),size(PFC.Classification_correct_chosen_color_each_belief,5)),...
    reshape(PFC.Classification_correct_chosen_color_each_belief(:,:,2,2,:),size(PFC.Classification_correct_chosen_color_each_belief,1),size(PFC.Classification_correct_chosen_color_each_belief,2),size(PFC.Classification_correct_chosen_color_each_belief,5)));

PFC.Other_belief_chosen_color=horzcat(reshape(PFC.Classification_correct_chosen_color_each_belief(:,:,1,2,:),size(PFC.Classification_correct_chosen_color_each_belief,1),size(PFC.Classification_correct_chosen_color_each_belief,2),size(PFC.Classification_correct_chosen_color_each_belief,5)),...
    reshape(PFC.Classification_correct_chosen_color_each_belief(:,:,2,1,:),size(PFC.Classification_correct_chosen_color_each_belief,1),size(PFC.Classification_correct_chosen_color_each_belief,2),size(PFC.Classification_correct_chosen_color_each_belief,5)));

%%
for t=1:size(LIP.Same_belief_chosen_color,1)
    p_LIP_same(t)=z_test_function_bootstrap(reshape(nanmean(LIP.Same_belief_chosen_color(t,:,:),2),size(LIP.Same_belief_chosen_color,3),1),1/2);
    p_LIP_other(t)=z_test_function_bootstrap(reshape(nanmean(LIP.Other_belief_chosen_color(t,:,:),2),size(LIP.Same_belief_chosen_color,3),1),1/2);
    this_dist=reshape(nanmean(LIP.Same_belief_chosen_color(t,:,:),2)-nanmean(LIP.Other_belief_chosen_color(t,:,:),2),size(LIP.Same_belief_chosen_color,3),1);
    p_LIP_same_other(t) = z_test_function_bootstrap(this_dist,0);
    
    p_FEF_same(t)=z_test_function_bootstrap(reshape(nanmean(FEF.Same_belief_chosen_color(t,:,:),2),size(FEF.Same_belief_chosen_color,3),1),1/2);
    p_FEF_other(t)=z_test_function_bootstrap(reshape(nanmean(FEF.Other_belief_chosen_color(t,:,:),2),size(FEF.Same_belief_chosen_color,3),1),1/2);
    this_dist=reshape(nanmean(FEF.Same_belief_chosen_color(t,:,:),2)-nanmean(FEF.Other_belief_chosen_color(t,:,:),2),size(LIP.Same_belief_chosen_color,3),1);
    p_FEF_same_other(t) = z_test_function_bootstrap(this_dist,0);
    
    p_PFC_same(t)=z_test_function_bootstrap(reshape(nanmean(PFC.Same_belief_chosen_color(t,:,:),2),size(PFC.Same_belief_chosen_color,3),1),1/2);
    p_PFC_other(t)=z_test_function_bootstrap(reshape(nanmean(PFC.Other_belief_chosen_color(t,:,:),2),size(PFC.Same_belief_chosen_color,3),1),1/2);
    this_dist=reshape(nanmean(PFC.Same_belief_chosen_color(t,:,:),2)-nanmean(PFC.Other_belief_chosen_color(t,:,:),2),size(LIP.Same_belief_chosen_color,3),1);
    p_PFC_same_other(t) = z_test_function_bootstrap(this_dist,0);
end


%%
figure

subplot(3,1,1)
hold on
shadedErrorBar(window_start_list,nanmean(nanmean(LIP.Same_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(LIP.Same_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(LIP.Same_belief_chosen_color(:,:,:),2),5,3)]',{'Color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list,nanmean(nanmean(LIP.Other_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(LIP.Other_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(LIP.Other_belief_chosen_color(:,:,:),2),5,3)]',{'--','Color',color_for_ROI(1,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_LIP_same_other(:),0.85,color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(window_start_list))])
plot_significance_level(window_start_list,p_LIP_other(:),0.9,color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(window_start_list))])
plot_significance_level(window_start_list,p_LIP_same(:),0.95,color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'left';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.1 1])
text(-575,0.85,'Within - Across')
text(-575,0.9,'Across')
text(-575,0.95,'Withins')
box off
xlabel('Time to response')
ylabel('Chosen color decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
subplot(3,1,2)
hold on

shadedErrorBar(window_start_list,nanmean(nanmean(FEF.Same_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(FEF.Same_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(FEF.Same_belief_chosen_color(:,:,:),2),5,3)]',{'Color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list,nanmean(nanmean(FEF.Other_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(FEF.Other_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(FEF.Other_belief_chosen_color(:,:,:),2),5,3)]',{'--','Color',color_for_ROI(2,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_FEF_same_other(:),0.85,color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(window_start_list))])
plot_significance_level(window_start_list,p_FEF_other(:),0.9,color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(window_start_list))])
plot_significance_level(window_start_list,p_FEF_same(:),0.95,color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'left';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.1 1])
text(-575,0.85,'Within - Across')
text(-575,0.9,'Across')
text(-575,0.95,'Withins')
box off
xlabel('Time to response')
ylabel('Chosen color decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])

subplot(3,1,3)
hold on

shadedErrorBar(window_start_list,nanmean(nanmean(PFC.Same_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(PFC.Same_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(PFC.Same_belief_chosen_color(:,:,:),2),5,3)]',{'Color',color_for_ROI(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list,nanmean(nanmean(PFC.Other_belief_chosen_color(:,:,:),2),3),[prctile(nanmean(PFC.Other_belief_chosen_color(:,:,:),2),95,3),prctile(nanmean(PFC.Other_belief_chosen_color(:,:,:),2),5,3)]',{'--','Color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_PFC_same_other(:),0.85,color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(window_start_list))])
plot_significance_level(window_start_list,p_PFC_other(:),0.9,color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(window_start_list))])
plot_significance_level(window_start_list,p_PFC_same(:),0.95,color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'left';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.1 1])
text(-575,0.85,'Within - Across')
text(-575,0.9,'Across')
text(-575,0.95,'Withins')
box off
xlabel('Time to response')
ylabel('Chosen color decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])


%%
figure

subplot(3,2,1)
hold on
shadedErrorBar(window_start_list,nanmean(nanmean(LIP.Classification_correct_chosen_color(:,:,:),2),3),[prctile(nanmean(LIP.Classification_correct_chosen_color(:,:,:),2),95,3),prctile(nanmean(LIP.Classification_correct_chosen_color(:,:,:),2),5,3)]',{'Color',color_for_ROI(1,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_chosen_color_LIP(:),0.95,color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.2 0.8])
box off
xlabel('Time to response')
ylabel('Chosen color decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])

subplot(3,2,2)
hold on
shadedErrorBar(window_start_list,nanmean(nanmean(LIP.Classification_belief_correct(:,:,:),2),3),[prctile(nanmean(LIP.Classification_belief_correct(:,:,:),2),95,3),prctile(nanmean(LIP.Classification_belief_correct(:,:,:),2),5,3)]',{'Color',color_for_ROI(1,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_belief_LIP(:),0.95,color_for_ROI(1,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.1 1])
box off
xlabel('Time to response')
ylabel('Expected template decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])

subplot(3,2,3)
hold on
shadedErrorBar(window_start_list,nanmean(nanmean(FEF.Classification_correct_chosen_color(:,:,:),2),3),[prctile(nanmean(FEF.Classification_correct_chosen_color(:,:,:),2),95,3),prctile(nanmean(FEF.Classification_correct_chosen_color(:,:,:),2),5,3)]',{'Color',color_for_ROI(2,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_chosen_color_FEF(:),0.95,color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.2 0.8])
box off
xlabel('Time to response')
ylabel('Chosen color decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])

subplot(3,2,4)
hold on
shadedErrorBar(window_start_list,nanmean(nanmean(FEF.Classification_belief_correct(:,:,:),2),3),[prctile(nanmean(FEF.Classification_belief_correct(:,:,:),2),95,3),prctile(nanmean(FEF.Classification_belief_correct(:,:,:),2),5,3)]',{'Color',color_for_ROI(2,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_belief_FEF(:),0.95,color_for_ROI(2,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.1 1])
box off
xlabel('Time to response')
ylabel('Expected template decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])


subplot(3,2,5)
hold on
shadedErrorBar(window_start_list,nanmean(nanmean(PFC.Classification_correct_chosen_color(:,:,:),2),3),[prctile(nanmean(PFC.Classification_correct_chosen_color(:,:,:),2),95,3),prctile(nanmean(PFC.Classification_correct_chosen_color(:,:,:),2),5,3)]',{'Color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_chosen_color_PFC(:),0.95,color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.2 0.8])
box off
xlabel('Time to response')
ylabel('Chosen color decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])

subplot(3,2,6)
hold on
shadedErrorBar(window_start_list,nanmean(nanmean(PFC.Classification_belief_correct(:,:,:),2),3),[prctile(nanmean(PFC.Classification_belief_correct(:,:,:),2),95,3),prctile(nanmean(PFC.Classification_belief_correct(:,:,:),2),5,3)]',{'Color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_belief_PFC(:),0.95,color_for_ROI(3,:),[0.05, 0.01, 0.05/(length(window_start_list))])

yl=yline(1/2,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'center';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xline(300,'-','Color',[0.25 0.75 0.25],'LineWidth',1)
xl=xline(200,'-',{'Response'},'Color',[0.25 0.75 0.25],'LineWidth',1);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
ylim([0.1 1])
box off
xlabel('Time to response')
ylabel('Expected template decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])



%%
function p = z_test_function_bootstrap(dist,null)

m = nanmean(dist);
s = nanstd(dist);
z = (m-null)/s;
p=1-normcdf(z);

end

function plot_significance_level(x,p,a,c,thr)
hold on
for b=1:length(thr)
    this_thr=thr(b);
    for i=1:length(p)-1
        if p(i)<=this_thr && p(i+1)<=this_thr
            plot(x(i):x(i+1),a*ones(1,x(i+1)-x(i)+1),'-','Color',c,'LineWidth',b)
        end
        if i>1
            if p(i-1)>=this_thr && p(i)<=this_thr && p(i+1)>=this_thr
                plot(x(i),a,'.','Color',c,'LineWidth',b)
            end
        end
        if p(end-1)>=this_thr && p(end)<=this_thr
            plot(x(end),a,'.','Color',c,'LineWidth',b)
        end
    end
end


end










