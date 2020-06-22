%% MRI ORIENTATION IN ACPC COORDSYS (non-interactive)
%% Virtual channel analysis of epilepsy MEG data
%% Kurtosis beamformer
function aTEST(a)
%%
restoredefaultpath %think to avoid issues always best to start from a tabula rasa
fieldtrip_path = '/home/c1725992/Downloads/fieldtrip-20200220/'; %don't need fieldtrip to be this up to date - can use latest CUBRIC version I believe
addpath(fieldtrip_path)
ft_defaults
addpath /cubric/software/MEG/CUBRIC/freqfun/ %for bandpassfilter.m (cheers Megan)
%% load in secondary files
participant_ID = importdata('/cubric/data/c1725992/BECTS_mrkd/aTEST_pipeline/participant_ID.txt'); % only difference is 2 datasets (091012-50 & 100211-51)
participant_dir = importdata('/cubric/data/c1725992/BECTS_mrkd/aTEST_pipeline/participant_dir.txt'); % naming convention can differ slightly so be careful
dataset_num = importdata('/cubric/data/c1725992/BECTS_mrkd/aTEST_pipeline/dataset_num.txt');
MEG_data = importdata('/cubric/data/c1725992/BECTS_mrkd/aTEST_pipeline/MEG_data.txt');
%% add looping variable to secondary files
participant_ID = participant_ID{a}; %there may be subtle differences in file naming convention
participant_dir = participant_dir{a}; %so it's best to be wary to avoid unnecessary errors
dataset_num = dataset_num{a};
MEG_data = MEG_data{a};
%%
datadir = ['/home/c1725992/data/BECTS_mrkd/' participant_dir]; %OG baseline data directory
%%
output_dir = ['/home/c1725992/data/BECTS_mrkd_prcssd/' participant_dir];
cd(output_dir);
MRI_data = fullfile(datadir, [participant_ID '_coreg.mri']); %OG baseline MRI file
mri = ft_read_mri(MRI_data); %can't change coordsys here, but can later on in script
%%
MEG_data = fullfile(datadir, MEG_data); %change the string each time - can use a looped variable that you create earlier in the script - see 'src_localisation.m'
%%
cfg = [];
ft_sourceplot(cfg, mri);
disp(mri);
close;
%%
headshape_file = fullfile(datadir, [participant_ID '_coreg.shape']); %need this information in order to change mri coordsys & orientation
headshape = ft_read_headshape(headshape_file);
disp(headshape);
%% non-interactive (so *fully* automatic)
cfg = [];
cfg.method = 'headshape';
cfg.headshape.headshape = headshape; 
cfg.headshape.interactive = 'no';
cfg.headshape.icp = 'no';
mri_realign = ft_volumerealign(cfg, mri);
ft_sourceplot([], mri_realign); 
disp(mri_realign);
%%
mri_realign = ft_convert_coordsys(mri_realign,'acpc'); %able to change the coordsys here but not sure why
disp(mri_realign);
%%
mri_resliced = ft_volumereslice([], mri_realign);
ft_sourceplot([], mri_resliced); 
disp(mri_resliced);

save mri_resliced.mat mri_resliced % save the data for sharing
%%
cfg = [];
cfg.dataset = MEG_data;
cfg.hpfilter = 'yes';
cfg.hpfreq = 20;
cfg.lpfilter = 'yes';
cfg.lpfreq = 70;
Bandpass = [20,70];
cfg.channel = 'MEG';
cfg.coilaccuracy = 0; % ensure that sensors are expressed in SI units
data = ft_preprocessing(cfg);

%Krish add
cfg=[];
cfg.detrend='no';
cfg.resamplefs=300; %Downsample to 300Hz for speed/memory
data=ft_resampledata(cfg,data);
%% neeed to change data coordsys
data.grad = ft_convert_coordsys(data.grad,'acpc');
disp(data.grad);
%%
cfg = [];
cfg.channel = 'MEG';
cfg.covariance = 'yes';
cov_matrix = ft_timelockanalysis(cfg, data);
disp(cov_matrix);
%% need to change cov_matrix coordsys
cov_matrix.grad=ft_convert_coordsys(cov_matrix.grad,'acpc');
disp(cov_matrix.grad); %does the coordsys look correct?
%%
if isfile ('seg.mat') == 1 % this variable takes a long time to be created
    load seg; % saves time loading in variable if you have it
    disp(seg);
elseif exist ('seg','var') == 0
cfg = [];
cfg.tissue = 'brain';
cfg.spmversion = 'spm12';
seg = ft_volumesegment(cfg, mri_resliced);
disp(seg);

save seg seg
end
%%
cfg = [];
cfg.tissue = 'brain';
cfg.spmversion = 'spm12';
brain_mesh = ft_prepare_mesh(cfg, seg);

cfg = [];
cfg.method = 'singleshell';
cfg.unit = 'm'; % ensure that the headmodel is expressed in SI units
headmodel = ft_prepare_headmodel(cfg, brain_mesh);
disp(headmodel);
%%
cfg = [];
cfg.method = 'basedonmri';
cfg.mri = mri_resliced;
cfg.headmodel = headmodel;
cfg.resolution = 0.004; % in SI units (7mm which is coarse) (5mm and below which is fine)
cfg.unit = 'm'; % ensure that the sourcemodel is expressed in SI units
cfg.grad = data.grad;
sourcemodel = ft_prepare_sourcemodel(cfg);
disp(sourcemodel);
%% visualize the source model, in combination with MRI & head model
% ensure all geometrical data properly aligned

figure
ft_plot_headmodel(headmodel, 'unit', 'mm');
ft_plot_sens(data.grad, 'unit', 'mm', 'coilsize', 10, 'chantype', 'meggrad');
ft_plot_mesh(sourcemodel.pos, 'unit', 'mm');
ft_plot_ortho(mri_resliced.anatomy, 'transform', mri_resliced.transform, 'style', 'intersect', 'unit', 'mm');
alpha 0.5

close;
%%
if isfile (['leadfield' dataset_num '.mat']) == 1
    load (['leadfield' dataset_num])
    disp(leadfield)
elseif exist ('leadfield','var') == 0
    cfg = [];
    cfg.channel = 'MEG';
    cfg.headmodel = headmodel;
    cfg.sourcemodel = sourcemodel;
    cfg.normalize = 'yes'; % normalization avoids power bias towards centre of head
    cfg.reducerank = 2;
    leadfield = ft_prepare_leadfield(cfg, cov_matrix);
    disp(leadfield);
    
save (['leadfield' dataset_num], 'leadfield')
end
%%
if isfile (['source' dataset_num '.mat']) == 1
    load (['source' dataset_num])
    disp(source)
elseif exist ('source','var') == 0
cfg = [];
cfg.headmodel = headmodel;
cfg.sourcemodel = leadfield;
cfg.method = 'lcmv';
cfg.lcmv.projectmom = 'yes'; % project dipole time series in direction of maximal power (see below)
cfg.lcmv.kurtosis = 'yes'; % compute kurtosis at each location
cfg.projectnoise  = 'yes';
cfg.lambda = '5%';
cfg.keepfilter    = 'yes'; % keeps the weights which you need to get virtual sensors​
source = ft_sourceanalysis(cfg, cov_matrix);
disp(source);

save (['source' dataset_num], 'source')
end
%%
% source is in m, mri_resliced is in mm, hence source_interp will also be in mm
cfg = [];
cfg.parameter = 'kurtosis';
source_interp = ft_sourceinterpolate(cfg, source, mri_resliced);
disp(source_interp);
%%
cfg = [];
cfg.filename = [participant_dir '_anatomy.nii'];
cfg.parameter = 'anatomy';
cfg.format = 'nifti';
ft_volumewrite(cfg, source_interp);

cfg = [];
cfg.filename = [participant_dir '_kurtosis' dataset_num '.nii'];
cfg.parameter = 'kurtosis';
cfg.format = 'nifti';
cfg.datatype = 'float'; % integer datatypes will be scaled to the maximum, floating point datatypes not
ft_volumewrite(cfg, source_interp);
%% Visualise the beamformer time series (want to use AnyWave software but don't know how to get it to work)
Fsample = data.fsample; % data variable is from ft_preprocessing
%% Get filters
filters = cell2mat(source.avg.filter(source.inside == 1)); %The filters come out backwards from fieldtrip for SAM vs LCMV! For LCMV no dash needed...
%% Normalise filters (vector norm as in Hillebrand 2012)
for i = 1 : size(filters,1)
        filters(i,:) = filters(i,:) ./ norm(filters(i,:));
end
%% Big concatanated loop 
NumTrials = length(data.trial);

LongSensor=[];

for t = 1:NumTrials
    LongSensor = [LongSensor double(data.trial{t})*1E15]; 
end
%% Get virtual sensors
if isfile (['VEs' dataset_num '.mat']) == 1
    load (['VEs' dataset_num]) 
    disp(VEs)
elseif exist ('VEs','var') == 0
VEs.matrix = bandpassfilter(filters*LongSensor, Fsample, Bandpass); % Bandpass is the frequency band you're using, e.g. [1,150] for 1-150Hz
VEs.fs = Fsample;
disp(VEs);

save (['VEs' dataset_num], 'VEs')
end
%% Sliding window (overlapping)
if isfile (['kurtcont' dataset_num '.mat']) == 1
    load (['kurtcont' dataset_num])
elseif exist ('kurtcont','var') == 0
sz=size(VEs.matrix, 2); %length of second dimension (columns/timepoints) of variable VE?
leng=2; %length of sliding window (in seconds) (adjustable)
jump=0.1; %length of incremental window jump (in seconds) (adjustable)
Ntim=VEs.fs*leng; %number of timepoints in window
Njum=VEs.fs*jump; %number of jumps?
Nwind=floor((sz-Ntim)/Njum); %number of windows
kurtcont=zeros(size(VEs.matrix, 1),Nwind); %matrix for the kurtosis values

for i=1:Nwind %from start to the end window
    slid1=1+(i-1)*Njum; %sliding window
    slid2 = slid1+Ntim;
    SOI=VEs.matrix(:,slid1:slid2); %every row value and the values between the sliding windows (signal of interest)
    givkurt=kurtosis(SOI,[],2)-3; %get kurtosis values ('-3' to make the kurtosis of the normal distribution equal to zero) 
    
  kurtcont(:,i)=givkurt; %save kurtosis values to matrix (in this context would check Nwind to estimate time of operation)
    if mod(i,250)==0
        
        fprintf(1,'Done %d windows out of %d \n',i,Nwind);

    end        
end

save (['kurtcont' dataset_num], 'kurtcont')
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         