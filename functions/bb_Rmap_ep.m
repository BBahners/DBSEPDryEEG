function [Rmap,Rperm,Rnaned,Rpermnaned,h, regressorperm]=bb_Rmap_ep(varargin)
% Create Correlation Matrix (=R-map or R-matrix)
% This function uses a Cell Array of the Evoked Potential (EP) time series
% (=fingerprints) consisting of channels x time matrices for each
% hemisphere or patient and loads each of the files to concatenate them
% into the Matrix X across all hemispheres or patients.Then this Matrix X
% is used to correlate it for every element of this matrix (element-wise
% correlation)with the individual patient's outcome (=regressor). The
% resulting correlation coefficients are put into a R-map (or R-matrix).
% This is based on the Lead DBS function ea_Rmap used for the DBS network 
% mapping approach by Horn et al. 2017. AnnNeurol.
% Author: Bahne H Bahners, BWH 2024

%usage:[Rmap,...]=bb_Rmap_ep(patConnectivityFiles(otherpts),...
% Regressor(otherpts),['Rmaps',filesep,'R_',sprintf('%02.0f',pt),'.nii'],...
% corrtype, 'save','permute');

fingerprints=varargin{1};
regressor=varargin{2};
output=varargin{3};
corrtype=varargin{4};
%saving=varargin{5};


pthresh=0.05;
itercount=5000;
h=nan;

dim1=size(fingerprints{1,1},1); % define dimension 1 (=channels)
dim2=size(fingerprints{1,1},2); % define dimension 2 (=time)
dim3=size(fingerprints,1); % define dimension 3 (= patients/hemispheres)
dim4=size(fingerprints{1,1}(:),1); % dimension4= long concatenated vector

if length(size(fingerprints{1,1}))>2 % in case of time-frequency map as input
dim5=size(fingerprints{1,1},3);
end

for ii=1:size(fingerprints,1)
    if ~exist('X','var')
        n=fingerprints{ii};
        X=n(:);
    else
        n=fingerprints{ii};
        X=[X,n(:)];
    end
end

X=X';
nnanix=~isnan(nansum(X,1)).*(abs(nansum(X,1)))>0;

R=corr(regressor,X,'type',corrtype,'rows','pairwise');

if exist('dim5')
Rmap=reshape(R,[dim1,dim2,dim5]);
else
Rmap=reshape(R,[dim1,dim2]);
end

if nargin>4
    if strcmp(varargin{5},'dontsave')
    else
        save(output,"Rmap");
    end

else
    save(output,"Rmap");
end

if nargin>5
    if ismember(varargin{6},{'permute','permuteplot'})
        Rperm=nan(itercount,size(X,2));
        regressorperm=repmat(regressor,1,itercount);
        for i=1:itercount
            regressorperm(:,i)=regressorperm(randperm(numel(regressor)),i);
        end
        Rperm(:,nnanix)=corr(regressorperm,X(:,nnanix),'type',corrtype,'rows','pairwise');
    end
end
if exist('Rperm','var') % permutation test
    Rpermnaned=[R;Rperm]; % for now, first entry is the unpermuted one.
    sRd=sort([R;Rperm],1,'descend');
    % delete values from Rpermnaned that are not significant (uncorrected):
    delp=Rpermnaned<...
        repmat(sRd(round((pthresh/2)*itercount),:),itercount+1,1);
    deln=Rpermnaned>...
        repmat(sRd(round((1-(pthresh/2))*itercount),:),itercount+1,1);
    del=logical(delp.*deln);
    Rpermnaned(del)=nan;

    Rnaned=Rpermnaned(1,:);
    Rpermnaned=Rpermnaned(2:end,:);




    if strcmp(varargin{6},'permuteplot')
        ea_dispercent(0,'Iterating patients');
        for pt=1:length(varargin{1})
            Xthispt=X(pt,:)';
            Ihat(pt)=atanh(corr(Xthispt,R','rows','pairwise','type',corrtype)); % real predictions
            Ihat_Rperm(pt,:)=atanh(corr(Xthispt,...
                Rperm','rows','pairwise','type',corrtype)); % permuted predictions
            ea_dispercent(pt/length(varargin{1}));
        end
        ea_dispercent(1,'end');

        [R]=corr(Ihat',regressor); % Predictive R of unpermuted values

        [Rperm]=corr(Ihat_Rperm,regressor); % Predictive R of permuted values

        Rperm=sort(Rperm,'descend'); % build distribution
        RlargerRperm=R>Rperm;
        p_predict_perm=sum(~RlargerRperm)./size(Rperm,1); % calculate final permutation based p value
        disp(['Permutation based p for overall prediction = ',num2str(p_predict_perm),'.']);

        h=ea_corrplot(regressor,Ihat',p_predict_perm,{'Empirical vs. Predicted','Empirical','Predicted'});
    end

end