%% Phantom Figures
clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
close all
% Machine-specific addpath removed; see config/local_paths.m


rmappath=dryeeg_result_path(repoRoot, 'phantom_plot_phantom_map');
if ~exist(rmappath,"file")
    mkdir(rmappath)
end

files={...
    dryeeg_data_path(repoRoot, 'PhantomComp/archive/P007_data_stim_average_240213_1817.mat');...
    dryeeg_data_path(repoRoot, 'PhantomComp/archive/P016_data_stim_average_240610_1455.mat');...
    dryeeg_data_path(repoRoot, 'PhantomComp/archive/Phantom_data_stim_average_240209_1759.mat')};
load(dryeeg_data_path(repoRoot, 'channels.mat')); 
channels(:,33:35)=[];
for ii=1:size(channels,2)
    channlabels{ii}=channels(ii).Name;
end
chanord=[sort(find(contains(string(char(channlabels)),'Fp')));find(contains(string(char(channlabels)),'AF'));find(contains(string(char(channlabels)),'F')&~contains(string(char(channlabels)),'A')&~contains(string(char(channlabels)),'Fp')&~contains(string(char(channlabels)),'FC'));find(contains(string(char(channlabels)),'FC'));find(contains(string(char(channlabels)),'C')&~contains(string(char(channlabels)),'F')&~contains(string(char(channlabels)),'P'));find(contains(string(char(channlabels)),'CP'));find(contains(string(char(channlabels)),'P')&~contains(string(char(channlabels)),'C'));find(contains(string(char(channlabels)),'T')&~contains(string(char(channlabels)),'P'));find(contains(string(char(channlabels)),'O')&~contains(string(char(channlabels)),'P'))];
chanlist_new=char(channlabels{chanord});
chanflip={'P7',	'T7','CP5',	'FC5',	'F7',	'F3',	'C3',	'P3',	'AF3',	'Fp1',	'Fp2',	'AF4',	'Fz',	'FC1',	'Cz',	'CP1',	'PO4',	'O2','Oz',	'O1',	'PO3',	'Pz',	'CP2',	'FC2',	'P4',	'C4',	'F4',	'F8',	'FC6',	'CP6',	'T8',	'P8'};
% Find the corresponding indices for channel reordering
[~, flip] = ismember(chanflip,channlabels);
for bb=1:length(files)
    orig= load(files{bb});
    if bb==1
        Patient=orig.F(1:32,:);
    elseif bb==2
        Patient2=orig.F(1:32,:);
    elseif bb==3
        Phantom=orig.F(1:32,:);
    end
end
EEGabs= (Phantom)*-10e6;
EMGabs=(Patient)*10e6;
MEGabs=(Patient2)*10e6;

cls=mandrill(3);
for pp=1:length(cls)
    c{pp}=cls(pp,:);
end
c{2}=[0.6,0.6,0.6];
time=orig.Time%;*1000;

h9=figure;
for ichanns=1:32
    mamp=MEGabs(ichanns,:);
    h(1)=plot(time,mamp,'Color',c{3},'LineWidth',2);
    hold on

end

for ichanns=1:32
    mamp=EMGabs(ichanns,:);
    h(2)=plot(time,mamp,'Color',c{1},'LineWidth',2);
    hold on
end

for ichanns=1:32
    if ichanns~=7
        mamp=EEGabs(ichanns,:);
        
        h(3)= plot(time,mamp,'Color',c{2},'LineWidth',2.5);
        hold on
       
    end
    hold on

    ax=gca;
    ax.XGrid='on';
    ax.YGrid='on';
    ax.XLabel.String='Time (s)';
    ax.XLim=[-0.02 0.15];
    ax.YLim=[-200 200];
    %ax.YScale='log';
    ax.YLabel.String='Amplitude (ÂµV)';
    ax.TickDir='out';
    ax.Box='off';
    ax.FontSize=10;
    fig=gcf;
    fig.Color='white';
end
legend(h,{'Patient 7','Patient 16','Phantom'})
ax.Legend.Box='off';
ax=gca;
ax.TickDir="out";
ax.Box="off";
ax.XGrid="on";
ax.YGrid="on";
ax.FontSize=14;
ax.FontWeight="bold";
saveas(h9,[rmappath,'PhantomvsPatients','.fig']);
saveas(h9,[rmappath,'PhantomvsPatients','.png']);

%% Figure 1
cmap=mandrill(256);
channelvec=1:32;

h2= figure
im=imagesc(time,channelvec,MEGabs(chanord,:))
colormap(cmap)
clim([-100,100]);
xlim([-0.03 0.15]);
xlabel('Time (s)','FontSize',10)
yticklabels(chanlist_new)
yticks(channelvec)
cb = colorbar();
ylabel(cb,'Amplitude (ÂµV)','FontSize',10,'Rotation',90)
cb.Label.Position(1) = 3;
cb.FontSize=14;
ax=gca;
ax.TickDir="out";
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight="bold";
ax.FontWeight='bold';
h2.Color='white';

saveas(h2,[rmappath,'PhantomvsPatientsColormap','.fig']);
 saveas(h2,[rmappath,'PhantomvsPatientsColormap','.png']);


MEGvec=(MEGabs);
EEGvec=(EEGabs);
mask=MEGvec>max(EEGvec)|MEGvec<min(EEGvec);
EEGabs_new=EEGabs;
EEGabs_new=EEGabs_new.^0;
EEGabs_new(~mask)=NaN;
h3= figure
im=imagesc(time,channelvec,EEGabs_new(chanord,:))
colormap('gray')
clim([0,1]);
xlim([-0.03 0.15]);
im.AlphaData=0.5;
xlabel('Time (s)','FontSize',10)
yticklabels(chanlist_new)
yticks(channelvec)
cb = colorbar();
ylabel(cb,'Amplitude (ÂµV)','FontSize',10,'Rotation',90)
cb.Label.Position(1) = 3;
cb.FontSize=14;
ax=gca;
ax.TickDir="out";
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight="bold";
ax.FontWeight='bold';
h3.Color='white';

 saveas(h3,[rmappath,'PhantomvsPatientsColormapMask','.fig']);
 saveas(h3,[rmappath,'PhantomvsPatientsColormapMask','.png']);

