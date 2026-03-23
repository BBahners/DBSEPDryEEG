%% Dry EEG Figures
% Author: Bahne Bahners
clear variables
addpath ./spm12/
addpath(genpath('./leaddbs'));
addpath /functions/
addpath /functions/plot_topography/ % Matlab function to plot topographies
 
load /functions/mandrillcolormap.mat
load /data/timevec.mat
load /data/channels.mat

allpts=1:29; % all patients vector (cave: only used for first loop), then hemispheres (allconds)
intp=1000; % interpolation for Topography plots

% Template plot for panel positions
Regressor=randn(1,29);
h= ea_corrplot(Regressor,Regressor);

%% Prepare Channel Order & Channel Flip
channels(:,33:35)=[];
for ii=1:size(channels,2)
    channlabels{ii}=channels(ii).Name; % create cell array with EEG channel labels
end
chanord= bb_chanord(channlabels);% now get a new channelorder for the plots (Frontal to Occipital)
chanlist_new=char(channlabels{chanord}); % now reorder the channellist
% new channellabels are assigned to essentially flip the EEG data
chanside={'right','right','right','right','right','right','right','right','right','right','left','left','z','right','z','right','left','left','z','right','right','z','left','left','left','left','left','left','left','left','left','left'};
chanflip={'P7',	'T7', 'CP5',	'FC5',	'F7',	'F3',	'C3',	'P3',	'AF3',	'Fp1',	'Fp2',	'AF4',	'Fz',	'FC1',	'Cz',	'CP1',	'PO4',	'O2',   'Oz',	'O1',	'PO3',	'Pz',	'CP2',	'FC2',	'P4',	'C4',	'F4',	'F8',	'FC6',	'CP6',	'T8',	'P8'};
[~, flip] = ismember(chanflip,channlabels); % Find the corresponding indices for channel reordering
ridx=contains(chanside,'right');
lidx=contains(chanside,'left');
channelvec=1:32;
allconds=1:58;


%% Figure 2
% Redefine EP map time and reload maps
% Analysis time window:
ini=-0.03; % beginning of time window in seconds
fin=0.2 ;% end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
time=time(win); % time based on window defined

load('/data/Amap.mat');
load('/data/Smap.mat');

% Time windows of interest
% Time window 1
init=0.020;
fint=0.030;
win1=[min(find(time>=init)):min(find(time>=fint))];
% Time window 2
init=0.05;
fint=0.06;
win2=[min(find(time>=init)):min(find(time>=fint))];
% Time window 3
init=0.075;
fint=0.085;
win3=[min(find(time>=init)):min(find(time>=fint))];
% Time window 4
init=0.090;
fint=0.12;
win4=[min(find(time>=init)):min(find(time>=fint))];

timwin={win1,win2,win3,win4};
timwinlab={'   20 to 30 ms', '   50 to 60 ms','   75 to 85 ms','   90 to 120 ms'};

% Figure 2 A
% Define Color maps
mand=mandrill(256);
mand=mand(end:-1:1,:);
channlabels_new=cellstr(regexprep(string(chanlist_new),' ',''));
vcmap=mandrill(size(Amap,1)+10);
vcmap=vcmap(end:-1:1,:);
vcmap(19:23,:)=[];

h4=figure;
for ch=1:size(Amap,1)
    plot(time,Amap(chanord(ch),:)*10e6,"Color",vcmap(ch,:))
    hold on
    conf= ciplot(Amap(chanord(ch),:)*10e6-Smap(chanord(ch),:)*10e6,Amap(chanord(ch),:)*10e6+Smap(chanord(ch),:)*10e6,time,vcmap(ch+5,:),0.2);
    conf.EdgeColor='none';
    ylim([-80,80]);
    xlim([ini,fin]);
end
xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (µV)','FontSize',14,'Rotation',90)
h4.Color='white';
ax=gca;
ax.TickDir="out";
ax.Box="off";
ax.XGrid="on";
ax.YGrid="on";
ax.FontSize=14;
ax.FontWeight="bold";
text(0.024, 67, 'P25','FontSize',14,'FontWeight','bold')
text(0.055, -40, 'N55','FontSize',14,'FontWeight','bold')
text(0.08, 15, 'P80','FontSize',14,'FontWeight','bold')
text(0.1, -20, 'N100','FontSize',14,'FontWeight','bold')

ax.Colormap=vcmap;
c=colorbar;

c.Ticks=[0:1/31:1];
c.Direction="reverse";
c.TickLabels=channlabels_new;
c.TickDirection="out";
h4.Position=h.Position;
h4.Position(4)=h4.Position(4)+100;

% Figure 2 B
% Topoplots
h3=figure;
for iwin=1:size(timwin,2)
    a= subplot(2,4,iwin);
    plot_topography(channlabels,(mean(Amap(1:32,timwin{iwin})')*10e6),false(1),'10-20',false(1),false(1),intp)
    fg=gcf;
    fg.Colormap=cmap;
    title(timwinlab{iwin},'FontSize',10);
    cb = colorbar();
    clim([-30 30])
    cb.Visible='off';
    ax=gca;
    ax.TickDir="out";
    ax.Children(1).LineWidth=1;
    ax.Children(2).LineWidth=1;
    ax.Children(3).LineWidth=1;
    ax.Children(4).LineWidth=1;
    a.Position(2)=a.Position(2)+0.1
    a.Position (1)=a.Position(1)-0.06
end

a=subplot(2,4,5:8); % A-Map

imagesc(time(min(find(time>=ini)):min(find(time>=fin))),[1:32],Amap(chanord,:)*10e6)
colormap(cmap)
clim([-30,30]);
xlabel('Time (s)','FontSize',10)
yticklabels(chanlist_new)
yticks([1:32])
cb = colorbar();
ylabel(cb,'Amplitude (µV)','FontSize',10,'Rotation',90)
cb.Label.Position(1) = 3;
cb.FontSize=14;
cb.FontWeight='bold';
h3.Color='white';
ax=gca;
ax.TickDir="out";
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight="bold";
ax.FontWeight='bold';

%Position move %% position =[x_position y_position widht length] all are in some unit
h3.Position=h.Position;
a.Position(4)=a.Position(4)+0.27 ;
a.Position(2)=a.Position(2)+0.015 ;
h3.Position(4)=h3.Position(4)+100


%% Figure 3
load('/data/Rmap.mat');
load('/data/timevec.mat');

%% Analysis time window:
ini=0.01; % beginning of time window in seconds
fin=0.2 ; % end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
Tmap=repmat(time(win),[32,1]); %time map
time=time(win); % time based on window defined

% Time windows of interest
% Time window 1
init=0.02;
fint=0.03;
win1=[min(find(time>=init)):min(find(time>=fint))];
% Time window 2
init=0.05;
fint=0.06;
win2=[min(find(time>=init)):min(find(time>=fint))];
% Time window 3
init=0.075;
fint=0.085;
win3=[min(find(time>=init)):min(find(time>=fint))];
% Time window 4
init=0.09;
fint=0.12;
win4=[min(find(time>=init)):min(find(time>=fint))];

timwin={win1,win2,win3,win4};
timwinlab={'   20 to 30 ms', '   50 to 60 ms','   75 to 85 ms','   90 to 120 ms'};


% Figure 3 B
% Topoplots
halfidx=contains(chanside,'left')|contains(chanside,'z') ;
h2=figure;
for iwin=1:size(timwin,2)

    a=subplot(2,4,iwin);

    [k,ch_x,ch_y]=plot_topography(channlabels,(mean(Rmap(1:32,timwin{iwin})')),false(1),'10-20',false(1),false(1),intp);%,'false')%,1,1,0)
    hold on
    scatter(ch_x(halfidx,:), ch_y(halfidx,:), 10,'white', 'LineWidth',1);

    fg=gcf;
    fg.Colormap=cmap;
    title(timwinlab{iwin},'FontSize',10);
    cb = colorbar();
    clim([-0.4 0.4]);
    cb.Visible='off';
    ax=gca;
    ax.TickDir="out";
    ax.Children(1).LineWidth=1;
    ax.Children(2).LineWidth=1;
    ax.Children(3).LineWidth=1;
    ax.Children(4).LineWidth=1;
    a.Position (2)=a.Position(2)+0.05;
    a.Position (1)=a.Position(1)-0.06;

end

a=subplot(2,4,5:8); % R-Matrix

Rmap(~halfidx,:)=NaN;
Rmap=Rmap(chanord,:);
chanlist_new(sum(isnan(Rmap),2)>1,:)=[];
channelvechalf=1:length(chanlist_new);
Rmap(sum(isnan(Rmap),2)>1,:)=[];

imagesc(time,channelvechalf,Rmap)
colormap(cmap)
clim([min(min(Rmap)),max(max(Rmap))]);
xlabel('Time (s)','FontSize',16)
yticklabels(chanlist_new)
yticks(channelvechalf)
cb = colorbar();
ylabel(cb,'R','FontSize',10,'Rotation',0)
cb.FontSize=14;
cb.FontWeight='bold';
h2.Color='white';
ax=gca;
ax.TickDir="out";
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight="bold";
ax.FontWeight='bold';

h2.Position=h.Position;
a.Position(4)=a.Position(4)+0.18;
a.Position(2)=a.Position(2)+0.03 ;




