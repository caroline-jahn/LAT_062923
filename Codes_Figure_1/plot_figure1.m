%plot cool stuff about the winning model
clear all

fsroot='/Volumes/buschman';
task='Learning_attentional_templates';

%set monkey
monkey='Beaker';

%do not change these settings
N_channels=6;
model_type='RW';
switch_type='Reset';
N_bins=100;

% load data
load_path = sprintf('/Volumes/buschman/Projects/Learning_Attentional_Templates/Analysed_data/%s',monkey);
load_name=sprintf('All_sessions_%s_%s_%d_channels_VBMC',switch_type,model_type,N_channels);
load(fullfile(load_path,load_name));
data=load_monkey_data_continuous(fsroot,monkey,N_channels,100);

PLOT_IND=1; %if you want to see all the blocks value functions (as in Fig 1G)

%plot option (do not change!)
N_plot=80;
N_plot_start=35;
load('colors')
color_current=colors(17,:); %green
color_old=colors(1,:); %pink

switch monkey
    case 'Beaker'
        cfp=[0.5 0.5 0.5];
    case 'Scooter'
        cfp=[0.5 0.5 0.5];
end

cd('/Volumes/buschman/Users/Caroline/Learning_attentional_templates/Behavioral_Analysis/Figure1')

%% Reorganise

block_sw_ind=find(Model_predictions.model_input.block==1);
day_sw_ind=find(data.X_data(:,5)==1 & data.X_data(:,7)==0);

for j=1:length(block_sw_ind)
    Best_chosen_for_plot(:,j)=Model_predictions.behavior.Best_chosen(block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    SecondBest_chosen_for_plot(:,j)=Model_predictions.behavior.SecondBest_chosen(block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    Worst_chosen_for_plot(:,j)=Model_predictions.behavior.Worst_chosen(block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    
    Best_chosen_predicted_for_plot(:,j)=Model_predictions.model_outputs.Choice_prediction(1,block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    SecondBest_chosen_predicted_for_plot(:,j)=Model_predictions.model_outputs.Choice_prediction(2,block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    Worst_chosen_predicted_for_plot(:,j)=Model_predictions.model_outputs.Choice_prediction(3,block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    
    RPE_for_plot(:,j)=Model_predictions.model_outputs.RPE(block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    Precision_for_plot(:,j)=Model_predictions.model_outputs.Precision(block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    Reward_for_plot(:,j)=data.Reward(block_sw_ind(j):block_sw_ind(j)+N_plot-1);
    
    Progression_for_plot(:,j)=data.X_data(block_sw_ind(j):block_sw_ind(j)+N_plot-1,5)./data.X_data(block_sw_ind(j):block_sw_ind(j)+N_plot-1,6);
end

for j=1:length(block_sw_ind)-1
    Switch_in_block(j)=sum(Model_predictions.model_outputs.Switch(block_sw_ind(j):block_sw_ind(j+1)-1));
end
Switch_in_block(length(block_sw_ind))=sum(Model_predictions.model_outputs.Switch(block_sw_ind(length(block_sw_ind)):end));

rule_angle=0:2*pi/100:2*pi-2*pi/100;

for j=1:length(block_sw_ind)
    new_rule(j)=data.X_data(block_sw_ind(j),1);
    if sum(day_sw_ind==block_sw_ind(j))==0
        old_rule(j)=data.X_data(block_sw_ind(j),8);
        for i=1:N_plot
            [~, Old_best(i,j)]=(min(abs(mod(data.X_data(i+block_sw_ind(j),2:4)-old_rule(j)+pi,2*pi)-pi)));
            Best_chosen_to_old(i,j)=data.choices(Old_best(i,j),i+block_sw_ind(j));
            Best_chosen_predicted_to_old(i,j)=Model_predictions.model_outputs.Choice_prediction(Old_best(i,j),i+block_sw_ind(j));
            
            Distance_rule_mean_belief_for_plot(i,j)=abs(mod(new_rule(j)-angle(sum(Model_predictions.model_outputs.Value_for_choice(:,i+block_sw_ind(j)).*exp(1i.*rule_angle')))+pi,2*pi)-pi);
        end
    else
        old_rule(j)=NaN;
        Old_best(1:N_plot,j)=NaN;
        Best_chosen_to_old(1:N_plot,j)=NaN;
        Best_chosen_predicted_to_old(1:N_plot,j)=NaN;
        for i=1:N_plot
            Distance_rule_mean_belief_for_plot(i,j)=abs(mod(new_rule(j)-angle(sum(Model_predictions.model_outputs.Value_for_choice(:,i+block_sw_ind(j)).*exp(1i.*rule_angle')))+pi,2*pi)-pi);
        end
    end
end

for j=1:length(block_sw_ind)
    if ~isnan(old_rule(j))
        for i=-N_plot:N_plot
            [~, best_ind]=min(abs(mod(new_rule(j)-data.options_color(:,i+block_sw_ind(j))+pi,2*pi)-pi));
            Min_error_to_best(i+N_plot+1,j)=abs(mod(new_rule(j)-data.options_color(best_ind,i+block_sw_ind(j))+pi,2*pi)-pi);
            Error_to_best(i+N_plot+1,j)=abs(mod(new_rule(j)-data.chosen_color(i+block_sw_ind(j))+pi,2*pi)-pi);
            Rel_error_to_best(i+N_plot+1,j)=Error_to_best(i+N_plot+1,j)-Min_error_to_best(i+N_plot+1,j);
            
            [~, best_ind]=min(abs(mod(old_rule(j)-data.options_color(:,i+block_sw_ind(j))+pi,2*pi)-pi));
            Min_error_to_old_best(i+N_plot+1,j)=abs(mod(old_rule(j)-data.options_color(best_ind,i+block_sw_ind(j))+pi,2*pi)-pi);
            Error_to_old_best(i+N_plot+1,j)=abs(mod(old_rule(j)-data.chosen_color(i+block_sw_ind(j))+pi,2*pi)-pi);
            Rel_error_to_old_best(i+N_plot+1,j)=Error_to_old_best(i+N_plot+1,j)-Min_error_to_old_best(i+N_plot+1,j);
        end
    else
        Min_error_to_best(1:2*N_plot+1,j)=NaN;
        Error_to_best(1:2*N_plot+1,j)=NaN;
        Rel_error_to_best(1:2*N_plot+1,j)=NaN;
        Min_error_to_old_best(1:2*N_plot+1,j)=NaN;
        Error_to_old_best(1:2*N_plot+1,j)=NaN;
        Rel_error_to_old_best(1:2*N_plot+1,j)=NaN;
    end
end

for j=1:length(block_sw_ind)
    if ~isnan(old_rule(j))
        Best_chosen_to_switch(:,j)=Model_predictions.behavior.Best_chosen(block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        SecondBest_chosen_to_switch(:,j)=Model_predictions.behavior.SecondBest_chosen(block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        Worst_chosen_to_switch(:,j)=Model_predictions.behavior.Worst_chosen(block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        
        Best_chosen_predicted_to_switch(:,j)=Model_predictions.model_outputs.Choice_prediction(1,block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        SecondBest_chosen_predicted_to_switch(:,j)=Model_predictions.model_outputs.Choice_prediction(2,block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        Worst_chosen_predicted_to_switch(:,j)=Model_predictions.model_outputs.Choice_prediction(3,block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        
        RPE_to_switch(:,j)=Model_predictions.model_outputs.RPE(block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        Precision_to_switch(:,j)=Model_predictions.model_outputs.Precision(block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
        for i=1:N_plot
            Distance_rule_mean_belief_to_switch(i,j)=abs(mod(old_rule(j)-angle(sum(Model_predictions.model_outputs.Value_for_choice(:,block_sw_ind(j)-N_plot+i-1).*exp(1i.*rule_angle')))+pi,2*pi)-pi);
        end
        Reward_to_switch(:,j)=data.Reward(block_sw_ind(j)-N_plot:block_sw_ind(j)-1);
    else
        Best_chosen_to_switch(1:N_plot,j)=NaN;
        SecondBest_chosen_to_switch(1:N_plot,j)=NaN;
        Worst_chosen_to_switch(1:N_plot,j)=NaN;
        
        Best_chosen_predicted_to_switch(1:N_plot,j)=NaN;
        SecondBest_chosen_predicted_to_switch(1:N_plot,j)=NaN;
        Worst_chosen_predicted_to_switch(1:N_plot,j)=NaN;
        
        RPE_to_switch(1:N_plot,j)=NaN;
        Precision_to_switch(1:N_plot,j)=NaN;
        
        Distance_rule_mean_belief_to_switch(1:N_plot,j)=NaN;
        Reward_to_switch(1:N_plot,j)=NaN;
    end
end

%%
for i=1:20
    p_bin(i)=myBinomTest(sum(Best_chosen_for_plot(i,:),2),size(Best_chosen_for_plot,2),1/3,'one');
end
Above_chance_choice=find(p_bin<0.05/20);
%%

rule_switch=abs(mod(new_rule-old_rule+pi,2*pi)-pi);
for j=1:length(block_sw_ind)
    
    Precision_start(j)=mean(Precision_for_plot(Progression_for_plot(:,j)<=1/3,j),1);
    Accuracy_start(j)=mean(Best_chosen_for_plot(1:N_plot_start,j),1);
    Accuracy_start_model(j)=mean(Best_chosen_predicted_for_plot(1:N_plot_start,j),1);
end

%% Figure S1 D

figure
subplot(1,2,2)
histogram(Model_predictions.model_outputs.Total_Switch_when,'FaceColor',cfp)
box off
title(sprintf('%s',monkey))
xlabel('Trials since switch','FontSize',20)
ylabel('Number of blocks with a reset at this trial','FontSize',20)
ax=gca;
ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;
xlim([0 50])
switch monkey
    case 'Beaker'
        set(gca,'YTick',(0:1:5));
end


%% Figure S1C

figure
subplot(1,2,2)
plot(rule_switch(Switch_in_block==2),Accuracy_start_model(Switch_in_block==2),'kd','MarkerFaceColor','k')
hold on
plot(rule_switch(Switch_in_block==1),Accuracy_start_model(Switch_in_block==1),'o','Color',[0.35 0.35 0.35],'MarkerFaceColor',[0.35 0.35 0.35])
hold on
plot(rule_switch(Switch_in_block==0),Accuracy_start_model(Switch_in_block==0),'o','Color',[0.75 0.75 0.75],'MarkerFaceColor',[0.75 0.75 0.75])
xlabel('Template switch (in rad)','FontSize',20)
ylabel(sprintf('Mean accuracy (model)\nin first %d trials',N_plot_start),'FontSize',20)
box off
% title(sprintf('%s',monkey))
set(gca,'XTick',(0:pi/4:pi));
set(gca,'XTickLabel',{'0','π/4','π/2','3π/2','π'})
ylim([0 1])
xlim([0 pi])
ax=gca;

ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

subplot(1,2,1)
plot(rule_switch(Switch_in_block==2),Accuracy_start(Switch_in_block==2),'kd','MarkerFaceColor','k')
hold on
plot(rule_switch(Switch_in_block==1),Accuracy_start(Switch_in_block==1),'o','Color',[0.35 0.35 0.35],'MarkerFaceColor',[0.35 0.35 0.35])
hold on
plot(rule_switch(Switch_in_block==0),Accuracy_start(Switch_in_block==0),'o','Color',[0.75 0.75 0.75],'MarkerFaceColor',[0.75 0.75 0.75])
xlabel('Template switch (in rad)','FontSize',20)
ylabel(sprintf('Mean accuracy (data)\nin first %d trials',N_plot_start),'FontSize',20)
box off
% title(sprintf('%s',monkey))
set(gca,'XTick',(0:pi/4:pi));
set(gca,'XTickLabel',{'0','π/4','π/2','3π/2','π'})
ylim([0 1])
xlim([0 pi])
ax=gca;

ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

%% correlation angular change and accuracy in first 35 trials

[r_m, pm]=corr(Accuracy_start_model(~isnan(rule_switch))',rule_switch(~isnan(rule_switch))');
[r_d, p_d]=corr(Accuracy_start(~isnan(rule_switch))',rule_switch(~isnan(rule_switch))');

%% Figure S1 E

no_reset=[0 0 0 0];
one_reset=[0 0 0 0];
two_reset=[0 0 0 0];

no_reset(1)=sum(Switch_in_block(rule_switch<=pi/4)==0);
one_reset(1)=sum(Switch_in_block(rule_switch<=pi/4)==1);
two_reset(1)=sum(Switch_in_block(rule_switch<=pi/4)==2);

no_reset(2)=sum(Switch_in_block(rule_switch>pi/4 & rule_switch<=pi/2)==0);
one_reset(2)=sum(Switch_in_block(rule_switch>pi/4 & rule_switch<=pi/2)==1);
two_reset(2)=sum(Switch_in_block(rule_switch>pi/4 & rule_switch<=pi/2)==2);

no_reset(3)=sum(Switch_in_block(rule_switch>pi/2 & rule_switch<=3*pi/4)==0);
one_reset(3)=sum(Switch_in_block(rule_switch>pi/2 & rule_switch<=3*pi/4)==1);
two_reset(3)=sum(Switch_in_block(rule_switch>pi/2 & rule_switch<=3*pi/4)==2);

no_reset(4)=sum(Switch_in_block(rule_switch>3*pi/4)==0);
one_reset(4)=sum(Switch_in_block(rule_switch>3*pi/4)==1);
two_reset(4)=sum(Switch_in_block(rule_switch>3*pi/4)==2);

figure
subplot(1,2,2)
bar(0.75:3.75,no_reset,0.25,'FaceColor',[0.75 0.75 0.75])
hold on
bar(1:4,one_reset,0.25,'FaceColor',[0.35 0.35 0.35])
hold on
bar(1.25:4.25,two_reset,0.25,'FaceColor','k')

xlabel(sprintf('Magnitude of \n template switch (in rad)'),'FontSize',16)
ylabel('Number of blocks')
box off
legend({'No reset', 'One reset', 'Two resets'},'FontSize',14,'Location','NorthWest')
% title(sprintf('%s',monkey))
set(gca,'XTick',(0.5:1:4.5));
set(gca,'XTickLabel',{'0','π/4','π/2','3π/2','π'})
ax=gca;

ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;


%% Figure 1 F

distance_mean_belief_rule_all=abs(mod(data.X_data(:,1)'-angle(sum(Model_predictions.model_outputs.Value_for_choice(:,:).*exp(1i.*rule_angle')))+pi,2*pi)-pi);

%Extract peak belief
N_bins=100;
color_binned =  0:(2*pi/N_bins):(2*pi-(2*pi/N_bins));

for i=1:size(Model_predictions.model_outputs.Value_for_choice,2)
    [~, peak_belief_index(i)]=max(Model_predictions.model_outputs.Value_for_choice(:,i));
end
peak_belief=color_binned(peak_belief_index);
distance_peak_belief_rule_all=abs(mod(data.X_data(:,1)'-peak_belief+pi,2*pi)-pi);

progression=data.X_data(:,5)./data.X_data(:,6);

for i=1:30
    distance_peak_belief_rule_all_binned(i)=mean(distance_peak_belief_rule_all(Model_predictions.model_outputs.Precision>=(i-1)/10 & Model_predictions.model_outputs.Precision<=(i)/10));
    distance_peak_belief_rule_all_binned_sem(i)=std(distance_peak_belief_rule_all(Model_predictions.model_outputs.Precision>=(i-1)/10 & Model_predictions.model_outputs.Precision<=(i)/10))/sqrt(sum(Model_predictions.model_outputs.Precision>=(i-1)/10 & Model_predictions.model_outputs.Precision<=(i)/10)-1);
end
for i=1:30
    distance_peak_belief_rule_all_binned_prog(i)=mean(distance_peak_belief_rule_all(progression>=(i-1)/30 & progression<=(i)/30));
    distance_peak_belief_rule_all_binned_prog_sem(i)=std(distance_peak_belief_rule_all(progression>=(i-1)/30 & progression<=(i)/30))/sqrt(sum(progression>=(i-1)/30 & progression<=(i)/30)-1);
end

figure
subplot(1,2,1)
shadedErrorBar([],distance_peak_belief_rule_all_binned_prog,distance_peak_belief_rule_all_binned_prog_sem,{'color',cfp,'LineWidth',2},2)
box off
ylabel(sprintf('Distance between template estimate \n and the true template (in rad)'),'FontSize',20)
xlabel('Progression in block','FontSize',20)
set(gca,'XTick',(0:10:30));
set(gca,'XTickLabel',{'0','1/3','2/3','1'})
ylim([0 pi/2+pi/8])
set(gca,'YTick',(0:pi/4:pi/2));
set(gca,'YTickLabel',{'0','π/4','π/2'})
ax=gca;

ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

Figure_1_F_mean=distance_peak_belief_rule_all_binned_prog;
Figure_1_F_sem=distance_peak_belief_rule_all_binned_prog_sem;


%% New vs. old rule
color_current_data=rgb2hsv(color_current);
color_current_data=hsv2rgb(color_current_data(1),color_current_data(2),0.65);

color_old_data=rgb2hsv(color_old);
color_old_data=hsv2rgb(color_old_data(1),color_old_data(2),0.9);

color_current_model=rgb2hsv(color_current);
color_current_model=hsv2rgb(color_current_model(1),color_current_model(2),0.35);

color_old_model=rgb2hsv(color_old);
color_old_model=hsv2rgb(color_old_model(1),color_old_model(2),0.5);

%% Figure 1 E

figure
subplot(1,2,2)
plot(nanmean(Best_chosen_for_plot,2),'Color',color_current_data,'Linewidth',2)
hold on
plot(nanmean(Best_chosen_to_old,2),'Color',color_old_data,'Linewidth',2)
hold on

plot(nanmean(Best_chosen_predicted_for_plot,2),'--','Linewidth',2,'Color',color_current_model)
hold on
plot(nanmean(Best_chosen_predicted_to_old,2),'--','Linewidth',2,'Color',color_old_model)
hold on

shadedErrorBar([],nanmean(Best_chosen_for_plot,2),nanstd(Best_chosen_for_plot,0,2)./sqrt(sum(~isnan(new_rule))-1),{'color',color_current_data},0.5)
hold on
shadedErrorBar([],nanmean(Best_chosen_to_old,2),nanstd(Best_chosen_to_old,0,2)./sqrt(sum(~isnan(old_rule))-1),{'color',color_old_data},0.5)
hold on
shadedErrorBar([],nanmean(Best_chosen_predicted_for_plot,2),nanstd(Best_chosen_predicted_for_plot,0,2)./sqrt(sum(~isnan(new_rule))-1),{'--','LineWidth',2,'color',color_current_model},0.5)
hold on
shadedErrorBar([],nanmean(Best_chosen_predicted_to_old,2),nanstd(Best_chosen_predicted_to_old,0,2)./sqrt(sum(~isnan(old_rule))-1),{'--','LineWidth',2,'color',color_old_model},0.5)

yline(0.33,'--')
legend({'Data: current template','Data: previous template', 'Model: current template','Model: previous template'})
ylim([0 1])
xlim([0 N_plot+1])
ylabel(sprintf('Proportion of trials in which\nthe best target was chosen'),'FontSize',20)
xlabel('Trial from template switch','FontSize',20)
box off
ax=gca;

ax.XAxis.FontSize=16;
ax.YAxis.FontSize=16;

Figure_1_E_current_data=Best_chosen_for_plot;
Figure_1_E_previous_data=Best_chosen_to_old;

Figure_1_E_current_model=Best_chosen_predicted_for_plot;
Figure_1_E_previous_model=Best_chosen_predicted_to_old;


%% Correlations model and behavior

[r_c, p_c]=corr(nanmean(Best_chosen_for_plot,2),nanmean(Best_chosen_predicted_for_plot,2));
[r_p, p_p]=corr(nanmean(Best_chosen_to_old,2),nanmean(Best_chosen_predicted_to_old,2));



%% Model correct predictions

mean_correct_prediction=mean(Model_predictions.model_outputs.Choice_prediction(data.choices==1));
std_correct_prediction=std(Model_predictions.model_outputs.Choice_prediction(data.choices==1));

%% Plot an example value function: figure 1G

if PLOT_IND==1
    
    for this_block=1:length(block_sw_ind)
        if ~isnan(old_rule(this_block))
            %find the bin of the current rule
            for c=1:N_bins
                if new_rule(this_block)>=(c-1)*2*pi/N_bins && new_rule(this_block)<=c*2*pi/N_bins
                    new_rule_binned=c;
                end
            end
            %find the bin of the old rule
            for c=1:N_bins
                if old_rule(this_block)>=(c-1)*2*pi/N_bins && old_rule(this_block)<=c*2*pi/N_bins
                    old_rule_binned=c;
                end
            end
            %find the bin of the expected template
            for j=-10:N_plot
                %peak belief
                if sum(Model_predictions.model_outputs.Value_for_choice(:,j+block_sw_ind(this_block))==0)==N_bins
                    PB(j+11)=NaN;
                else
                    PB(j+11)=mod(peak_belief(j+block_sw_ind(this_block)),2*pi);
                end
            end
            for i=1:length(PB)
                for c=1:N_bins
                    if ~isnan(PB(i))
                        if PB(i)>=(c-1)*2*pi/N_bins && PB(i)<=c*2*pi/N_bins
                            PB_binned(i)=c;
                        end
                    else
                        PB_binned(i)=NaN;
                    end
                end
            end
            
            %find the bin of the chosen color
            for j=-10:N_plot
                %chosen color
                CC(j+11)=mod(data.chosen_color(j+block_sw_ind(this_block)),2*pi);
            end
            for i=1:length(CC)
                for c=1:N_bins
                    if ~isnan(CC(i))
                        if CC(i)>=(c-1)*2*pi/N_bins && CC(i)<=c*2*pi/N_bins
                            CC_binned(i)=c;
                        end
                    else
                        CC_binned(i)=NaN;
                    end
                end
            end
            
            %find the bin of the other colors
            for j=-10:N_plot
                %other colors
                if data.Chosen_ind(j+block_sw_ind(this_block))==1
                    oc_id=[2,3];
                elseif data.Chosen_ind(j+block_sw_ind(this_block))==2
                    oc_id=[1,3];
                else
                    oc_id=[1,2];
                end
                OC(j+11,1)=mod(data.options_color(oc_id(1),j+block_sw_ind(this_block)),2*pi);
                OC(j+11,2)=mod(data.options_color(oc_id(2),j+block_sw_ind(this_block)),2*pi);
            end
            for i=1:length(CC)
                for c=1:N_bins
                    if OC(i,1)>=(c-1)*2*pi/N_bins && OC(i,1)<=c*2*pi/N_bins
                        OC_binned(i,1)=c;
                    end
                    if OC(i,2)>=(c-1)*2*pi/N_bins && OC(i,2)<=c*2*pi/N_bins
                        OC_binned(i,2)=c;
                    end
                end
            end
            this_value=zeros(size(Model_predictions.model_outputs.Value_for_choice,1),N_plot);
            for j=-10:N_plot
                for i=1:size(Model_predictions.model_outputs.Value_for_choice,1)
                    this_value(i,j+11)=Model_predictions.model_outputs.Value_for_choice(i,j+block_sw_ind(this_block))-min(Model_predictions.model_outputs.Value_for_choice(:,[-10:N_plot]+block_sw_ind(this_block)),[],'all');
                end
            end
            this_value=this_value./max(this_value(:,end));
            
            %color
            angle_rule=0:2*pi/N_bins:2*pi-2*pi/N_bins;
            for i=1:size(Model_predictions.model_outputs.Value_for_choice,1)
                for c=1:40
                    if angle_rule(i)>=(c-1)*2*pi/40 && angle_rule(i)<=(c)*2*pi/40
                        color_binned(i)=c;
                    end
                end
            end
            
            f=figure('visible','off');
            %f=figure('visible','on');
            for j=-10:N_plot
                for i=1:size(Model_predictions.model_outputs.Value_for_choice,1)
                    
                    x=[j-1+11 j+11 j+11 j-1+11];
                    y=[i-1 i-1 i i];
                    this_color=rgb2hsv(colors(color_binned(i),:));
                    this_val=this_value(i,j+11)*this_color(3);
                    this_val(this_val>1)=1;
                    this_color=hsv2rgb([this_color(1),this_color(2),this_val]);
                    patch(x,y,this_color,'EdgeColor',this_color)
                    hold on
                end
            end
            for j=-10:N_plot
                if Model_predictions.model_outputs.Switch(j+block_sw_ind(this_block))==1
                    xline(j+11-1,'r','LineWidth',2)
                    hold on
                end
                x=[j-1+11 j+11 j+11 j-1+11];
                %plot the peak
                if ~isnan(PB_binned(j+11))
                    y=[PB_binned(j+11)-1 PB_binned(j+11)-1 PB_binned(j+11) PB_binned(j+11)];
                    patch(x,y,this_color,'EdgeColor','w','FaceColor','w','LineWidth',2)
                    hold on
                end
                
                %plot the chosen color
                if ~isnan(CC_binned(j+11))
                    y=[CC_binned(j+11)-1 CC_binned(j+11)-1 CC_binned(j+11) CC_binned(j+11)];
                    patch(x,y,this_color,'EdgeColor','w','FaceColor','k','LineWidth',2)
                    hold on
                end
                
                %plot the other colors
                if ~isnan(OC_binned(j+11,1))
                    for k=1:2
                    y=[OC_binned(j+11,k)-1 OC_binned(j+11,k)-1 OC_binned(j+11,k) OC_binned(j+11,k)];
                    patch(x,y,this_color,'EdgeColor',[0.5 0.5 0.5],'FaceColor','k','LineWidth',2)
                    hold on
                    end
                end
                
            end
            hold on
            yline(new_rule_binned-1.5,'--','Color',colors(color_binned(new_rule_binned),:),'LineWidth',7)
            hold on
            yline(old_rule_binned-1.5,'--','Color',colors(color_binned(old_rule_binned),:),'LineWidth',7)
            
            set(gca,'YTick',(0:N_bins/4:N_bins));
            set(gca,'YTickLabel',{'0','π/2','π','3π/3','2π'})
            xlim([0 N_plot+11])
            set(gca,'XTick',(1:10:N_plot+11));
            set(gca,'XTickLabel',{'-10','0','10','20','30','40','50','60','70','80'})
            ax=gca;
            
            ax.XAxis.FontSize=16;
            ax.YAxis.FontSize=16;
            ylabel('Color','Fontsize',20)
            xlabel('Trials since switch','FontSize',20)
            
            clear *binned PB CC OC
            
        else
            this_value=zeros(size(Model_predictions.model_outputs.Value_for_choice,1),N_plot);
            %find the bin of the current rule
            for c=1:N_bins
                if new_rule(this_block)>=(c-1)*2*pi/N_bins && new_rule(this_block)<=c*2*pi/N_bins
                    new_rule_binned=c;
                end
            end
            %find the bin of the estimated template
            for j=1:N_plot
                if sum(Model_predictions.model_outputs.Value_for_choice(:,j+block_sw_ind(this_block))==0)==N_bins
                    PB(j)=NaN;
                else
                    PB(j)=mod(peak_belief(j+block_sw_ind(this_block)),2*pi);
                end
            end
            for i=1:length(PB)
                for c=1:N_bins
                    if ~isnan(PB(i))
                        if PB(i)>=(c-1)*2*pi/N_bins && PB(i)<=c*2*pi/N_bins
                            PB_binned(i)=c;
                        end
                    else
                        PB_binned(i)=NaN;
                    end
                end
            end
            %find the bin of the chosen color
            for j=1:N_plot
                %chosen color
                CC(j)=mod(data.chosen_color(j+block_sw_ind(this_block)),2*pi);
            end
            for i=1:length(CC)
                for c=1:N_bins
                    if ~isnan(CC(i))
                        if CC(i)>=(c-1)*2*pi/N_bins && CC(i)<=c*2*pi/N_bins
                            CC_binned(i)=c;
                        end
                    else
                        CC_binned(i)=NaN;
                    end
                end
            end
            %find the bin of the other colors
            for j=1:N_plot
                %other colors
                if data.Chosen_ind(j+block_sw_ind(this_block))==1
                    oc_id=[2,3];
                elseif data.Chosen_ind(j+block_sw_ind(this_block))==2
                    oc_id=[1,3];
                else
                    oc_id=[1,2];
                end
                OC(j,1)=mod(data.options_color(oc_id(1),j+block_sw_ind(this_block)),2*pi);
                OC(j,2)=mod(data.options_color(oc_id(2),j+block_sw_ind(this_block)),2*pi);
            end
            for i=1:length(CC)
                for c=1:N_bins
                    if OC(i,1)>=(c-1)*2*pi/N_bins && OC(i,1)<=c*2*pi/N_bins
                        OC_binned(i,1)=c;
                    end
                    if OC(i,2)>=(c-1)*2*pi/N_bins && OC(i,2)<=c*2*pi/N_bins
                        OC_binned(i,2)=c;
                    end
                end
            end
            for j=1:N_plot
                for i=1:size(Model_predictions.model_outputs.Value_for_choice,1)
                    this_value(i,j)=Model_predictions.model_outputs.Value_for_choice(i,j+block_sw_ind(this_block))-min(Model_predictions.model_outputs.Value_for_choice(:,[1:N_plot]+block_sw_ind(this_block)),[],'all');
                end
            end
            this_value=this_value./max(this_value(:,end));
            
            %color
            angle_rule=0:2*pi/N_bins:2*pi-2*pi/N_bins;
            for i=1:size(Model_predictions.model_outputs.Value_for_choice,1)
                for c=1:40
                    if angle_rule(i)>=(c-1)*2*pi/40 && angle_rule(i)<=(c)*2*pi/40
                        color_binned(i)=c;
                    end
                end
            end
            
            f=figure('visible','off');
            %f=figure('visible','on');
            for j=1:N_plot
                for i=1:size(Model_predictions.model_outputs.Value_for_choice,1)
                    
                    x=[j-1 j j j-1];
                    y=[i-1 i-1 i i];
                    this_color=rgb2hsv(colors(color_binned(i),:));
                    this_val=this_value(i,j)*this_color(3);
                    this_val(this_val>1)=1;
                    this_color=hsv2rgb([this_color(1),this_color(2),this_val]);
                    patch(x,y,this_color,'EdgeColor',this_color)
                    hold on
                end
            end
            
            for j=1:N_plot
                if Model_predictions.model_outputs.Switch(j+block_sw_ind(this_block))==1
                    xline(j-1,'r','LineWidth',2)
                    hold on
                end
                x=[j-1 j j j-1];
                if ~isnan(PB_binned(j))
                    y=[PB_binned(j)-1 PB_binned(j)-1 PB_binned(j) PB_binned(j)];
                    patch(x,y,this_color,'EdgeColor','w','FaceColor','w','LineWidth',2)
                    hold on
                end
                %plot the chosen color
                if ~isnan(CC_binned(j))
                    y=[CC_binned(j)-1 CC_binned(j)-1 CC_binned(j) CC_binned(j)];
                    patch(x,y,this_color,'EdgeColor','w','FaceColor','k','LineWidth',2)
                    hold on
                end
                
                %plot the other colors
                if ~isnan(OC_binned(j,1))
                    for k=1:2
                    y=[OC_binned(j,k)-1 OC_binned(j,k)-1 OC_binned(j,k) OC_binned(j,k)];
                    patch(x,y,this_color,'EdgeColor',[0.5 0.5 0.5],'FaceColor','k','LineWidth',2)
                    hold on
                    end
                end
            end
            
            hold on
            yline(new_rule_binned-1.5,'--','Color',colors(color_binned(new_rule_binned),:),'LineWidth',7)
            xlim([0 N_plot])
            set(gca,'YTick',(0:N_bins/4:N_bins));
            set(gca,'YTickLabel',{'0','π/2','π','3π/3','2π'})
            set(gca,'XTick',(1:10:N_plot));
            set(gca,'XTickLabel',{'0','10','20','30','40','50','60','70','80'})
            ax=gca;
            
            ax.XAxis.FontSize=16;
            ax.YAxis.FontSize=16;
            ylabel('Color (rad)','FontSize',20)
            xlabel('Trials since session start','FontSize',20)
            
            clear *binned PB CC OC
            
        end
        
        saveas(f,fullfile('/Volumes/buschman/Projects/Learning_Attentional_Templates/Analysed_data',sprintf('%s',monkey),'figures',sprintf('%d_peak_dashed_CC',this_block)),'epsc')
        
    end
    
end


