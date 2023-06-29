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

window_size=200;subsubtask_classifier=sprintf('Classifier_%s_%d_to_%s_%d',event,window_start_list(1),event,window_start_list(end));

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

N_boot=100;

dirstem = fullfile('Projects',task,'Analysis','Electrophy_analysis',subtask,subsubtask_classifier);
data_path_clasifier = fullfile(fsroot,dirstem);

color_for_ROI=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125]; %LIP FEF PFC

load('colors')

color_belief(1,:)=colors(1,:);
color_belief(2,:)=colors(14,:);
color_belief(3,:)=colors(27,:);

N_loc=1;

N_prog=1;

save_name='Chosen_value_peak_template_classifier_results_plot';

LOAD=0;


%%
if LOAD==1
    %LIP
    ROI='LIP';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            results_name=sprintf('Pseudo_pop_value_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);
            
            load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_value_each_belief_correct','Classification_value_correct','W_value','W','W_value_each_belief')
            
            for k=1:N_prog
                LIP.Classification_correct_belief(this_time,:,k,n_tt)=Classification_belief_correct(:,k);
                
                LIP.Classification_correct_value(this_time,:,k,n_tt)=Classification_value_correct(:,k);
                
                LIP.Classification_correct_value_each_belief(this_time,:,:,:,k,n_tt)=Classification_value_each_belief_correct(:,:,:,k);
            end
            LIP.W_template(this_time,:,:,n_tt)=W(:,:);
            LIP.W_value(this_time,:,n_tt)=W_value(:);
            LIP.W_value_each_belief(this_time,:,:,n_tt)=W_value_each_belief(:,:);
            
            clear Classification_belief_correct Classification_value_correct Classification_value_each_belief_correct W*
            
        end
    end
    
    %FEF
    ROI='FEF';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            results_name=sprintf('Pseudo_pop_value_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);
            
            load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_value_each_belief_correct','Classification_value_correct','W_value','W','W_value_each_belief')
            
            for k=1:N_prog
                FEF.Classification_correct_belief(this_time,:,k,n_tt)=Classification_belief_correct(:,k);
                
                FEF.Classification_correct_value(this_time,:,k,n_tt)=Classification_value_correct(:,k);
                
                FEF.Classification_correct_value_each_belief(this_time,:,:,:,k,n_tt)=Classification_value_each_belief_correct(:,:,:,k);
            end
            FEF.W_template(this_time,:,:,n_tt)=W(:,:);
            FEF.W_value(this_time,:,n_tt)=W_value(:);
            FEF.W_value_each_belief(this_time,:,:,n_tt)=W_value_each_belief(:,:);
            
            clear Classification_belief_correct Classification_value_correct Classification_value_each_belief_correct W*
            
        end
    end
    
    %PFC
    ROI='PFC';
    
    for n_tt=1:N_boot
        
        for this_time=1:length(window_start_list)
            
            results_name=sprintf('Pseudo_pop_value_peak_template_classifier_results_%s_%d_%d_%d',ROI,this_time,n_tt);
            
            load(fullfile(data_path_clasifier,results_name),'Classification_belief_correct','Classification_value_each_belief_correct','Classification_value_correct','W_value','W','W_value_each_belief')
            
            for k=1:N_prog
                PFC.Classification_correct_belief(this_time,:,k,n_tt)=Classification_belief_correct(:,k);
                
                PFC.Classification_correct_value(this_time,:,k,n_tt)=Classification_value_correct(:,k);
                
                PFC.Classification_correct_value_each_belief(this_time,:,:,:,k,n_tt)=Classification_value_each_belief_correct(:,:,:,k);
            end
            PFC.W_template(this_time,:,:,n_tt)=W(:,:);
            PFC.W_value(this_time,:,n_tt)=W_value(:);
            PFC.W_value_each_belief(this_time,:,:,n_tt)=W_value_each_belief(:,:);
            
            clear Classification_belief_correct Classification_value_correct Classification_value_each_belief_correct W*
            
        end
    end
    
    
    
    %%
    for i=1:length(window_start_list)
        for k=1:N_prog
            p_belief_LIP(i,k)=z_test_function_bootstrap(mean(LIP.Classification_correct_belief(i,:,k,:),2),1/3);
            p_belief_FEF(i,k)=z_test_function_bootstrap(mean(FEF.Classification_correct_belief(i,:,k,:),2),1/3);
            p_belief_PFC(i,k)=z_test_function_bootstrap(mean(PFC.Classification_correct_belief(i,:,k,:),2),1/3);
            
            p_value_LIP(i,k)=z_test_function_bootstrap(mean(LIP.Classification_correct_value(i,:,k,:),2),1/2);
            p_value_FEF(i,k)=z_test_function_bootstrap(mean(FEF.Classification_correct_value(i,:,k,:),2),1/2);
            p_value_PFC(i,k)=z_test_function_bootstrap(mean(PFC.Classification_correct_value(i,:,k,:),2),1/2);
            for l=1:3
                for m=1:3
                    p_value_each_belief_LIP(i,l,m,k)=z_test_function_bootstrap(mean(LIP.Classification_correct_value_each_belief(i,:,l,m,k,:),2),1/2);
                    p_value_each_belief_FEF(i,l,m,k)=z_test_function_bootstrap(mean(FEF.Classification_correct_value_each_belief(i,:,l,m,k,:),2),1/2);
                    p_value_each_belief_PFC(i,l,m,k)=z_test_function_bootstrap(mean(PFC.Classification_correct_value_each_belief(i,:,l,m,k,:),2),1/2);
                end
            end
        end
    end
    %%
    LIP.Same_belief_value=horzcat(reshape(LIP.Classification_correct_value_each_belief(:,:,1,1,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)),...
        reshape(LIP.Classification_correct_value_each_belief(:,:,2,2,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)),...
        reshape(LIP.Classification_correct_value_each_belief(:,:,3,3,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)));
    
    LIP.Other_belief_value=horzcat(reshape(LIP.Classification_correct_value_each_belief(:,:,1,2,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)),...
        reshape(LIP.Classification_correct_value_each_belief(:,:,1,3,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)),...
        reshape(LIP.Classification_correct_value_each_belief(:,:,2,1,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)),...
        reshape(LIP.Classification_correct_value_each_belief(:,:,2,3,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)),...
        reshape(LIP.Classification_correct_value_each_belief(:,:,3,1,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)),...
        reshape(LIP.Classification_correct_value_each_belief(:,:,3,2,1,:,:),size(LIP.Classification_correct_value_each_belief,1),size(LIP.Classification_correct_value_each_belief,2),size(LIP.Classification_correct_value_each_belief,6),size(LIP.Classification_correct_value_each_belief,7)));
    
    FEF.Same_belief_value=horzcat(reshape(FEF.Classification_correct_value_each_belief(:,:,1,1,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)),...
        reshape(FEF.Classification_correct_value_each_belief(:,:,2,2,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)),...
        reshape(FEF.Classification_correct_value_each_belief(:,:,3,3,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)));
    
    FEF.Other_belief_value=horzcat(reshape(FEF.Classification_correct_value_each_belief(:,:,1,2,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)),...
        reshape(FEF.Classification_correct_value_each_belief(:,:,1,3,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)),...
        reshape(FEF.Classification_correct_value_each_belief(:,:,2,1,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)),...
        reshape(FEF.Classification_correct_value_each_belief(:,:,2,3,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)),...
        reshape(FEF.Classification_correct_value_each_belief(:,:,3,1,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)),...
        reshape(FEF.Classification_correct_value_each_belief(:,:,3,2,1,:,:),size(FEF.Classification_correct_value_each_belief,1),size(FEF.Classification_correct_value_each_belief,2),size(FEF.Classification_correct_value_each_belief,6),size(FEF.Classification_correct_value_each_belief,7)));
    
    PFC.Same_belief_value=horzcat(reshape(PFC.Classification_correct_value_each_belief(:,:,1,1,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)),...
        reshape(PFC.Classification_correct_value_each_belief(:,:,2,2,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)),...
        reshape(PFC.Classification_correct_value_each_belief(:,:,3,3,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)));
    
    PFC.Other_belief_value=horzcat(reshape(PFC.Classification_correct_value_each_belief(:,:,1,2,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)),...
        reshape(PFC.Classification_correct_value_each_belief(:,:,1,3,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)),...
        reshape(PFC.Classification_correct_value_each_belief(:,:,2,1,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)),...
        reshape(PFC.Classification_correct_value_each_belief(:,:,2,3,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)),...
        reshape(PFC.Classification_correct_value_each_belief(:,:,3,1,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)),...
        reshape(PFC.Classification_correct_value_each_belief(:,:,3,2,1,:,:),size(PFC.Classification_correct_value_each_belief,1),size(PFC.Classification_correct_value_each_belief,2),size(PFC.Classification_correct_value_each_belief,6),size(PFC.Classification_correct_value_each_belief,7)));
    
    %%
    for t=1:size(LIP.Same_belief_value,1)
        p_LIP_same(t)=z_test_function_bootstrap(reshape(mean(LIP.Same_belief_value(t,:,:),2),size(LIP.Same_belief_value,3),1),1/2);
        p_LIP_other(t)=z_test_function_bootstrap(reshape(mean(LIP.Other_belief_value(t,:,:),2),size(LIP.Same_belief_value,3),1),1/2);
        this_dist(:)=mean(LIP.Same_belief_value(t,:,:),2)-mean(LIP.Other_belief_value(t,:,:),2);
        p_LIP_same_other(t) = z_test_function_bootstrap(this_dist,0);
        
        p_FEF_same(t)=z_test_function_bootstrap(reshape(mean(FEF.Same_belief_value(t,:,:),2),size(FEF.Same_belief_value,3),1),1/2);
        p_FEF_other(t)=z_test_function_bootstrap(reshape(mean(FEF.Other_belief_value(t,:,:),2),size(FEF.Same_belief_value,3),1),1/2);
        this_dist(:)=mean(FEF.Same_belief_value(t,:,:),2)-mean(FEF.Other_belief_value(t,:,:),2);
        p_FEF_same_other(t) = z_test_function_bootstrap(this_dist,0);
        
        p_PFC_same(t)=z_test_function_bootstrap(reshape(mean(PFC.Same_belief_value(t,:,:),2),size(PFC.Same_belief_value,3),1),1/2);
        p_PFC_other(t)=z_test_function_bootstrap(reshape(mean(PFC.Other_belief_value(t,:,:),2),size(PFC.Same_belief_value,3),1),1/2);
        this_dist(:)=mean(PFC.Same_belief_value(t,:,:),2)-mean(PFC.Other_belief_value(t,:,:),2);
        p_PFC_same_other(t) = z_test_function_bootstrap(this_dist,0);
    end
    
    save(fullfile(data_path_clasifier,save_name))
    
else
    
    load(fullfile(data_path_clasifier,save_name))
    
end


%%
for i=1:3
    color_for_ROI(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI(i,:)=hsv2rgb(color_for_ROI(i,1),color_for_ROI(i,2),1);
    color_for_ROI_late(i,:)=rgb2hsv(color_for_ROI(i,:));
    color_for_ROI_late(i,:)=hsv2rgb(color_for_ROI_late(i,1),color_for_ROI_late(i,2),0.55);
    
end


%%
figure
subplot(3,1,1)
hold on

plot(window_start_list,mean(mean(LIP.Classification_correct_value(:,:,1,:),2),4),'Color',color_for_ROI(1,:),'LineWidth',1)

shadedErrorBar(window_start_list,mean(mean(LIP.Classification_correct_value(:,:,1,:),2),4),[prctile(mean(LIP.Classification_correct_value(:,:,1,:),2),95,4),prctile(mean(LIP.Classification_correct_value(:,:,1,:),2),5,4)]',{'Color',color_for_ROI(1,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_value_LIP(:),1,color_for_ROI(1,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
ylim([0.35 1.05])
xl=xline(300,'-',{'Reward'},'Color',[0.25 0 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xl=xline(400,'-','Color',[0.25 0 0.75],'LineWidth',1);
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])



subplot(3,1,2)
hold on
plot(window_start_list,mean(mean(FEF.Classification_correct_value(:,:,1,:),2),4),'Color',color_for_ROI(2,:),'LineWidth',1)

shadedErrorBar(window_start_list,mean(mean(FEF.Classification_correct_value(:,:,1,:),2),4),[prctile(mean(FEF.Classification_correct_value(:,:,1,:),2),95,4),prctile(mean(FEF.Classification_correct_value(:,:,1,:),2),5,4)]',{'Color',color_for_ROI(2,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_value_FEF(:),1,color_for_ROI(2,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
ylim([0.35 1.05])
xl=xline(300,'-',{'Reward'},'Color',[0.25 0 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xl=xline(400,'-','Color',[0.25 0 0.75],'LineWidth',1);
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])


subplot(3,1,3)
hold on
plot(window_start_list,mean(mean(PFC.Classification_correct_value(:,:,1,:),2),4),'Color',color_for_ROI(3,:),'LineWidth',1)

shadedErrorBar(window_start_list,mean(mean(PFC.Classification_correct_value(:,:,1,:),2),4),[prctile(mean(PFC.Classification_correct_value(:,:,1,:),2),95,4),prctile(mean(PFC.Classification_correct_value(:,:,1,:),2),5,4)]',{'Color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_value_PFC(:),1,color_for_ROI(3,:),[0.01,0.05/(length(window_start_list)), 0.01/(length(window_start_list))])
yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
ylim([0.35 1.05])
xl=xline(300,'-',{'Reward'},'Color',[0.25 0 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xl=xline(400,'-','Color',[0.25 0 0.75],'LineWidth',1);
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])



%%
figure
subplot(3,1,1)
hold on


shadedErrorBar(window_start_list,mean(mean(LIP.Same_belief_value(:,:,:),2),3),[prctile(mean(LIP.Same_belief_value(:,:,:),2),95,3),prctile(mean(LIP.Same_belief_value(:,:,:),2),5,3)]',{'Color',color_for_ROI(1,:),'LineWidth',2},2)
shadedErrorBar(window_start_list,mean(mean(LIP.Other_belief_value(:,:,:),2),3),[prctile(mean(LIP.Other_belief_value(:,:,:),2),95,3),prctile(mean(LIP.Other_belief_value(:,:,:),2),5,3)]',{'--','Color',color_for_ROI(1,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_LIP_same_other(:),0.9,color_for_ROI(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list,p_LIP_other(:),0.95,color_for_ROI(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list,p_LIP_same(:),1,color_for_ROI(1,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
ylim([0.35 1.05])
xl=xline(300,'-',{'Reward'},'Color',[0.25 0 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xl=xline(400,'-','Color',[0.25 0 0.75],'LineWidth',1);
text(-575,0.9,'Within - across')
text(-575,0.95,'Across')
text(-575,1,'Within')
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])
subplot(3,1,2)
hold on
shadedErrorBar(window_start_list,mean(mean(FEF.Same_belief_value(:,:,:),2),3),[prctile(mean(FEF.Same_belief_value(:,:,:),2),95,3),prctile(mean(FEF.Same_belief_value(:,:,:),2),5,3)]',{'Color',color_for_ROI(2,:),'LineWidth',2},2)
shadedErrorBar(window_start_list,mean(mean(FEF.Other_belief_value(:,:,:),2),3),[prctile(mean(FEF.Other_belief_value(:,:,:),2),95,3),prctile(mean(FEF.Other_belief_value(:,:,:),2),5,3)]',{'--','Color',color_for_ROI(2,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_FEF_same_other(:),0.9,color_for_ROI(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list,p_FEF_other(:),0.95,color_for_ROI(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list,p_FEF_same(:),1,color_for_ROI(2,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
ylim([0.35 1.05])
xl=xline(300,'-',{'Reward'},'Color',[0.25 0 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xl=xline(400,'-','Color',[0.25 0 0.75],'LineWidth',1);
text(-575,0.9,'Within - across')
text(-575,0.95,'Across')
text(-575,1,'Within')
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])


subplot(3,1,3)
hold on
shadedErrorBar(window_start_list,mean(mean(PFC.Same_belief_value(:,:,:),2),3),[prctile(mean(PFC.Same_belief_value(:,:,:),2),95,3),prctile(mean(PFC.Same_belief_value(:,:,:),2),5,3)]',{'Color',color_for_ROI(3,:),'LineWidth',2},2)
shadedErrorBar(window_start_list,mean(mean(PFC.Other_belief_value(:,:,:),2),3),[prctile(mean(PFC.Other_belief_value(:,:,:),2),95,3),prctile(mean(PFC.Other_belief_value(:,:,:),2),5,3)]',{'--','Color',color_for_ROI(3,:),'LineWidth',2},2)

plot_significance_level(window_start_list,p_PFC_same_other(:),0.9,color_for_ROI(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list,p_PFC_other(:),0.95,color_for_ROI(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])
plot_significance_level(window_start_list,p_PFC_same(:),1,color_for_ROI(3,:),[0.01,0.05/length(window_start_list), 0.01/length(window_start_list)])

yl=yline(0.5,'--',{'Chance level'},'Color',[0 0 0],'LineWidth',1);
xline(-300,'-','Color',[0.5 0.5 0.5],'LineWidth',1)
xl=xline(-400,'-',{'Fixation'},'Color',[0.5 0.5 0.5],'LineWidth',1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'right';
yl.FontSize=12;
ylim([0.35 1.05])
xl=xline(300,'-',{'Reward'},'Color',[0.25 0 0.75],'LineWidth',1);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'left';
xl.FontSize=12;
xl=xline(400,'-','Color',[0.25 0 0.75],'LineWidth',1);
text(-575,0.9,'Within - across')
text(-575,0.95,'Across')
text(-575,1,'Within')
box off
xlabel('Time to targets on')
ylabel('Chosen value decoding accuracy')
xlim([-600 window_start_list(end)+150])
xticks([-600:200:window_start_list(end)+150])

%%

function p = z_test_function_bootstrap(dist,null)

m = mean(dist);
s = std(dist);
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
    end
    if p(end-1)>=this_thr && p(end)<=this_thr
        plot(x(end),a,'.','Color',c,'LineWidth',b)
    end
end


end










