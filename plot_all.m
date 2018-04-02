% Plots the steady-state lariat levels (defined as lariat reads mapped per
% million hg19 reads) for all datasets.  Includes patient fibroblast RNAseq
% data from all used genotypes (Controls, DBR1, TLR3, STAT1) with the all 
% treatments (no treatment, ifna2b stimulation, pic stimulation, 24hrs HSV1
% infection).  Plots each RNAseq experiment as a data point (each 
% experiment is an average of 3 technical replicates) and calculates error
% bars using SEM.
% (Lariat reads mapped using my lariat read aligner:
% https://github.com/allisontaggart/lariat, hg19 reads mapped using STAR w/
% default parameters, and counting done with in-house scripts)

% Full data/analysis as part of :
% "Inborn Errors of RNA Lariat Metabolism in Humans with Brainstem Viral
% Infection" - Cell. 2018 Feb 22;172(5):952-965.e18. doi: 10.1016/j.cell.2018.02.019.


close all; clear all;

set(0,'DefaultAxesFontSize',14);
set(0,'DefaultAxesFontName', 'Calibri')

%all_data.txt contains lariat levels (lariat reads mapped per million hg19
%reads) for all datasets including:
%   Control patients : no infection, no stimulation
%   Patient 1 : no infection, no stimulation
%   Patient 5 : no infection, no stimulation
%   Patient 6 : no infection, no stimulation
%   TLR3-/- Patient : no infection, no stimulation
%   STAT1-/- Patient : no infection, no stimulation
%   Control patients : ifna2b stimulation
%   Patient 1 : ifna2b stimulation
%   Patient 5 : ifna2b stimulation
%   Patient 6 : ifna2b stimulation
%   TLR3-/- Patient : ifna2b stimulation
%   STAT1-/- Patient : ifna2b stimulation
%   Control patients : pIC stimulation
%   Patient 1 : pIC stimulation
%   Patient 5 : pIC stimulation
%   Patient 6 : pIC stimulation
%   TLR3-/- Patient : pIC stimulation
%   STAT1-/- Patient : pIC stimulation
%   Control patients : 24h HSV1 infection
%   Patient 1 : 24h HSV1 infection
%   Patient 5 : 24h HSV1 infection
%   Patient 6 : 24h HSV1 infection
%   TLR3-/- Patient : 24h HSV1 infection
%   STAT1-/- Patient : 24h HSV1 infection

data = dlmread('./demo_data/all_data.txt','\t');
all = data(:,1);

%setting tick mark positions
x1 = [0.7 1.0 1.3];
x2 = [3 3.0 3];
x3 = [5 5.0 5];
x4 = [7 7.0 7];
x5 = [8.7 9.0 9.5];
x6 = [10.7 11.0 11.5];

% read in data (no stimuation, no infection)
ctl1 = all(4:6,:);
ctl2 = all(10:12,:);
ctl3 = all(16:18,:);
p1 = all(22:24,:);
p5 = all(28:30,:);
p6 = all(34:36,:);
s=all(40:42,:);
t = all(46:48,:);

% read in data (ifna2b stimulation)
ctl1_i = all(49:51,:);
ctl2_i = all(55:57,:);
ctl3_i = all(61:63,:);
p1_i = all(67:69,:);
p5_i = all(73:75,:);
p6_i = all(79:81,:);
s_i=all(85:87,:);
t_i = all(91:93,:);

% read in data (polyIC stimulation)
ctl1_p = all(52:54,:);
ctl2_p = all(58:60,:);
ctl3_p = all(64:66,:);
p1_p = all(70:72,:);
p5_p = all(76:78,:);
p6_p = all(82:84,:);
s_p=all(88:90,:);
t_p = all(94:96,:);

% read in data (24h HSV1 infection)
ctl1_24 = all(1:3,:);
ctl2_24 = all(7:9,:);
ctl3_24 = all(13:15,:);
p1_24 = all(19:21,:);
p5_24 = all(25:27,:);
p6_24 = all(31:33,:);
s_24 = all(37:39,:);
t_24 = all(43:45,:);


%set figure parameters
hFig = figure(1);
set(gcf,'Units','Inches');
set(hFig, 'Position', [.25 2.5 15 8])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Subplot - No immune stimulation/infection

subplot(2,2,1)

title('No Immune Stimulation/Infection');

hold on;

xlim([0 13]);
ylabel({'lariat reads / million';'mapped hg19 reads'})

%plot all data points with correct shape/color for genetic background
%plot data points individually, and then add SEM error bars

plot(x1,ctl1,'o','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl2,'o','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl3,'o','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',4,'MarkerEdgeColor',[0.5 0.5 0.5]);
ctl_sem=std(vertcat(ctl1,ctl2,ctl3))/sqrt(9);
errorbar(1,mean(vertcat(ctl3,ctl2,ctl1)),ctl_sem,'k');

plot(x6,s,'ko','MarkerFaceColor','k','MarkerSize',4);
s_sem = std(s)/sqrt(3);
errorbar(11,mean(s),s_sem,'k');

plot(x5,t,'ko','MarkerFaceColor','k','MarkerSize',4);
t_sem=std(t)/sqrt(3);
errorbar(9,mean(t),t_sem,'k');

plot(x2,p1,'ro','MarkerFaceColor', 'r','MarkerSize',4);
p1_sem=std(p1)/sqrt(3);
errorbar(3,mean(p1),p1_sem,'k');

plot(x3,p5,'o','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p5_sem=std(p5)/sqrt(3);
errorbar(5,mean(p5),p5_sem,'k');

plot(x4,p6,'^','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p6_sem=std(p6)/sqrt(3);
errorbar(7,mean(p6),p6_sem,'k');

%set labels
ax=gca;
ax.XTick= [1, 3, 5, 7, 9, 11];
ax.XTickLabels={'C' 'P1', 'P5', 'P6','TLR3-/-','STAT1-/-'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Subplot - IFN-a2b Stimulation
subplot(2,2,2)

title('IFN-a2b Stimulation');
hold on;
xlim([0 13]);


ylabel({'lariat reads / million';'mapped hg19 reads'})

%plot all data points with correct shape/color for genetic background
%plot data points individually, and then add SEM error bars
plot(x1,ctl1_i,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl2_i,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl3_i,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
ctl_i_sem=std(vertcat(ctl1_i,ctl2_i,ctl3_i)/sqrt(9));
errorbar(1,mean(vertcat(ctl1_i,ctl2_i,ctl3_i)),ctl_i_sem,'k');

plot(x6,s_i,'ko','MarkerFaceColor','k','MarkerSize',4);
s_i_sem = std(s_i)/sqrt(3);
errorbar(11,mean(s_i),s_i_sem,'k');

plot(x5,t_i,'ko','MarkerFaceColor','k','MarkerSize',4);
t_i_sem = std(t_i)/sqrt(3);
errorbar(9,mean(t_i),t_i_sem,'k');

plot(x2,p1_i,'ro','MarkerFaceColor', 'r','MarkerSize',4);
p1_i_sem=std(p1_i)/sqrt(3);
errorbar(3,mean(p1_i),p1_i_sem,'k');

plot(x3,p5_i,'o','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p5_i_sem=std(p5_i)/sqrt(3);
errorbar(5,mean(p5_i),p5_i_sem,'k');

plot(x4,p6_i,'^','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p6_i_sem=std(p6_i)/sqrt(3);
errorbar(7,mean(p6_i),4,p6_i_sem,'k');

%set labels
ax=gca;
ax.XTick= [1, 3, 5, 7, 9, 11];
ax.XTickLabels={'C' 'P1', 'P5', 'P6','TLR3-/-','STAT1-/-'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Subplot - polyIC Stimulation
subplot(2,2,3)

title('poly(I:C) Stimulation');

hold on;
xlim([0 13]);

ylabel({'lariat reads / million';'mapped hg19 reads'})

%plot all data points with correct shape/color for genetic background
%plot data points individually, and then add SEM error bars
plot(x1,ctl1_p,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl2_p,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl3_p,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
ctl_p_sem=std(vertcat(ctl1_p,ctl2_p,ctl3_p)/sqrt(9));
errorbar(1,mean(vertcat(ctl1_p,ctl2_p,ctl3_p)),ctl_p_sem,'k');

plot(x6,s_p,'ko','MarkerFaceColor','k','MarkerSize',4);
s_p_sem = std(s_p)/sqrt(3);
errorbar(11,mean(s_p),s_p_sem,'k');

plot(x5,t_p,'ko','MarkerFaceColor','k','MarkerSize',4);
t_p_sem = std(t_p)/sqrt(3);
errorbar(9,mean(t_p),t_p_sem,'k');

plot(x2,p1_p,'ro','MarkerFaceColor', 'r','MarkerSize',4);
p1_p_sem=std(p1_p)/sqrt(3);
errorbar(3,mean(p1_p),p1_p_sem,'k');

plot(x3,p5_p,'o','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p5_p_sem=std(p5_p)/sqrt(3);
errorbar(5,mean(p5_p),p5_p_sem,'k');

plot(x4,p6_p,'^','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p6_p_sem=std(p6_p)/sqrt(3);
errorbar(7,mean(p6_p),4,p6_p_sem,'k');

%set labels
ax=gca;
ax.XTick= [1, 3, 5, 7, 9, 11];
ax.XTickLabels={'C' 'P1', 'P5', 'P6','TLR3-/-','STAT1-/-'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Subplot - 24 hrs HSV1 infection

subplot(2,2,4)

title('24 hrs HSV1 Infection');
hold on;

xlim([0 13]);
ylim([0 150]);

ylabel({'lariat reads / million';'mapped hg19 reads'})

%plot all data points with correct shape/color for genetic background
%plot data points individually, and then add SEM error bars
plot(x1,ctl1_24,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl2_24,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
plot(x1,ctl3_24,'ko','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize', 4,'MarkerEdgeColor',[0.5 0.5 0.5]);
ctl_24_sem=std(vertcat(ctl1_24,ctl2_24,ctl3_24)/sqrt(9));
errorbar(1,mean(vertcat(ctl1_24,ctl2_24,ctl3_24)),ctl_24_sem,'k');

plot(x6,s_24,'ko','MarkerFaceColor','k','MarkerSize',4);
s_24_sem = std(s_24)/sqrt(3);
errorbar(11,mean(s_24),s_24_sem,'k');

plot(x5,t_24,'ko','MarkerFaceColor','k','MarkerSize',4);
t_24_sem = std(t_24)/sqrt(3);
errorbar(9,mean(t_24),t_24_sem,'k');

plot(x2,p1_24,'ro','MarkerFaceColor', 'r','MarkerSize',4);
p1_24_sem=std(p1_24)/sqrt(3);
errorbar(3,mean(p1_24),p1_24_sem,'k');

plot(x3,p5_24,'o','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p5_24_sem=std(p5_24)/sqrt(3);
errorbar(5,mean(p5_24),p5_24_sem,'k');

plot(x4,p6_24,'^','color',[0 0.5 0],'MarkerFaceColor', [0 0.5 0],'MarkerSize',4);
p6_24_sem=std(p6_24)/sqrt(3);
errorbar(7,mean(p6_24),4,p6_24_sem,'k');

%set labels
ax=gca;
ax.XTick= [1, 3, 5, 7, 9, 11];
ax.XTickLabels={'C' 'P1', 'P5', 'P6','TLR3-/-','STAT1-/-'};


