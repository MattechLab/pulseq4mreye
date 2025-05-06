%% This file intergrate Nils' phyllotaxis calculation
% The way to calculate timing will be different
% adapted by Yiwei Jia
% LIBRE GRE RADIAL PHYLLOTAXIS
%info:
% TR: 85.80ms 3.9ms fov 160mm(should double) TE 1.79ms flip angle 12
% nSeg:22 nShot:3275 matrix_size: 80 voxel_size:2 Slice_per_slab 80 (should be double)
% BW: 906Hz/px
% https://github.com/QIS-MRI/Pulseq/blob/main/write_pulseq_PhaseCycledbSSFP

rmpath(genpath('/Users/mauroleidi/Desktop/MattechGit/pulseSeqYiwei/pulseq'));rmpath(genpath('/Users/mauroleidi/Desktop/MattechGit/pulseSeqYiwei/pulseq'));
%%
clc;close all;clear;addpath(genpath('/Users/mauroleidi/Desktop/MattechGit/pulseSeqYiwei/pulseq'))
addpath(genpath('/Users/mauroleidi/Desktop/MattechGit/pulseSeqYiwei/pulseq4mreye/func'))
addpath(genpath('/Users/mauroleidi/Desktop/MattechGit/pulseSeqYiwei/pulseq/matlab'))
%% mode control
subject_num_array = 4;
subject_num = subject_num_array;
grad_mode = 'Whisper';
checking = 0;
overlapping_k_center = 1;

if subject_num ==3
    % gre_traj
    rf_spoiling_increment = 50; 
    adc_offset = 0;
    traj_design = 0;    % Trajectory design: 0 = original, 1 = pole-to-pole, 2 = continuous， 3=unishuffle
    notesuffix = 'oriTraj';
elseif subject_num ==4
    % gre_unishuffle
    rf_spoiling_increment = 50; 
    adc_offset = 0;
    traj_design = 3;    % Trajectory design: 0 = original, 1 = pole-to-pole, 2 = continuous, 3=unishuffle
    notesuffix = 'unishuTraj';
end


seq_mode = 2; %1: pre 2: main 3: debug
add_rfdelay = true;
rfdelay_ratio = 1;
do_compile = 1;
seqfolder = "/Users/mauroleidi/Desktop/MattechGit/pulseSeqYiwei/pulseq4mreye/mreye2p0_0430/";
use_rfspoiler = 1;

%%%%
seq_mode_list = {'pre', 'main', 'debug'};
seq_note = strcat('sub',num2str(subject_num),'_mauro_gre','_',seq_mode_list{seq_mode},'_', notesuffix);
disp(['Preparing: ', seq_note])
%%
flagSelfNav = 1;
nSeg = 22;
nShot_plot = 2;

if seq_mode == 1
    nShot = 419;
    nLine = nSeg*nShot;
    iLine = nLine;
    disp('--------------Prescan mode----------------')
elseif seq_mode == 2
    nShot = 3275;
    nLine = nSeg*nShot;
    iLine = nLine;
    disp('--------------Main sequence mode----------------')
else
    nShot = 102;
    nLine = nSeg*nShot;
    iLine = nLine;
     disp('--------------Debug sequence mode----------------')
end

%MatrixSize = 240*2; 
% The FoV is a parameter that depends on the size of the imaged region, to
% avoid the distruption of the image by artifacts we often double the FoV.
% So the AcquisionFoV = 2* ROIFoV
% Then it follows delta_k = 1/AcquisionFoV
% Finally we have to decide the value of N points per readouts, sometimes
% called matrix size in Siemens implementations.
% N sets the maximum sampled frequency in the k space = delta_k*N/2 =
% = N/(2*AcquisionFoV) This value is also known as Receiver bandwidth of the
% acquisition, and if too low can cause gibbs ringing artifacts

% Single FoV = 160 mm => AcquisionFoV = 320 mm
fov= [320*1e-3 320*1e-3 320*1e-3];     % Define FOV, isotropic, unit[m]
deltak=1 ./ fov;

% Points per line = To understand what is the correct number here
N = 240*2;
Nx      = N; %1e3*fov(1)/(res(2)); 
Ny      = N; %1e3*fov(2)/(res(2)); 
Nz      = N; %1e3*fov(3)/res(3);  
bw_px = 906; %Hz/px

%if checking
%    disp(['fov: ', num2str(fov)])
%end
% 1/bw: 
%res = [2 2 2];                         Define resolution, voxel size (mm)

%
switch grad_mode
    case 'Fast'
        max_grad = 26;      % Max gradient strength [mT/m]
        max_slew = 180;  % Maximum slew rate [mT/m/ms] orginial 120
    case 'Normal'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 100;     % Maximum slew rate [mT/m/ms]
    case 'Whisper'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 50;      % Maximum slew rate [mT/m/ms]
    case 'Performance'
        max_grad = 37;      % Max gradient strength [mT/m]
        max_slew = 188;  % Maximum slew rate [mT/m/ms]
end

sys = mr.opts('MaxGrad',max_grad,'GradUnit','mT/m',...
    'MaxSlew',max_slew,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6, 'B0', 2.89);

seq=mr.Sequence(sys);           % Create a new sequence object

%% Define the parameters

adc_dur = 1e6/bw_px; % Here adc duration seems fixed by the bandwidth, however it means that the bandwidth is per readout not per k_space sample
adc_dwell = round(adc_dur/N/(sys.rfRasterTime*1e6))*(sys.rfRasterTime*1e6); % Here the dwell time is basically the closest approximation of adc_dur/N, but a multiple of the rfRasterTime
adc_dur = adc_dwell*N; % Here we recompute the adc_duration but this time is a multiple of the rfRasterTime

%% WE SHOULD MAYBE TRY TO UNDERSTAND THE BANDWIDTH A BIT BETTER. Especially since I changed N, 
% and it impacts the adc_dwell I would like to understand what is the
% practical effect of that on the sequence, and if I should adapt other
% parameters accordingly
% Probably the sequence parameter should be the single sample bandwith,
% instead of the full readout bandwith? 



%% new gradient from Nils
roDuration = ceil((adc_dur*1e-6)/seq.gradRasterTime)*seq.gradRasterTime;
gz               = mr.makeTrapezoid('z','FlatArea',(Nx-1)*deltak(1),'FlatTime',roDuration,'system',sys);
adc  = mr.makeAdc(Nx,'Duration',gz.flatTime,'Delay',gz.riseTime, 'system',sys);

if overlapping_k_center
    gzPre_area = -gz.area/2-1/2*deltak(1);
else
    gzPre_area = -gz.area/2;
end
gzPre            = mr.makeTrapezoid('z','Area',gzPre_area,'system',sys);
prephase_min_dur = mr.calcDuration(gzPre);

gz_parts          = mr.splitGradientAt(gz,ceil(mr.calcDuration(adc)/sys.gradRasterTime)*sys.gradRasterTime);
%gz_parts{1}: rise+flat
%gz_parts{2}: down
%
gz_parts(1).delay = mr.calcDuration(gzPre);
%gz_1: gz_pre + gz rise +gz flat
gz_1              = mr.addGradients({gzPre,gz_parts(1)},'system',sys);

adc.delay         = adc.delay+mr.calcDuration(gzPre); 
% so from here adc starts with gz_1 (gz_pre)

%gz_2: gz down + gz spoiler
gz_parts(2).delay = 0;

% here different from nils'
gzSpoil=mr.makeTrapezoid('z','Area',0.5*Nx*deltak(3),'system',sys);
gzSpoil.delay       = mr.calcDuration(gz_parts(2));
gz_2              = mr.addGradients({gz_parts(2),gzSpoil},'system',sys);

gzSpoil.delay       = 0; 
pe_dur            = mr.calcDuration(gz_2);


%% RF objects with rf spoiler
% Creating the RF object
% How to determine the duration of rf pulse?
alpha = 12;                               % Define flip angle
df = -500;
tao = 1100; %us


% Spoiling increment in degrees
rf_spoiling_phase = 0;
adc_phaseOffset = [];
rf_spoiling_phase_rad = [];

for indLine = 1:iLine
    rf_spoiling_phase = indLine*(indLine - 1)*(rf_spoiling_increment/2);
    rf_spoiling_phase_rad{indLine} = rf_spoiling_phase * pi / 180;
    adc_phaseOffset{indLine} = rf_spoiling_phase_rad{indLine} + adc_offset;
    if mod(indLine, 500)==0
        disp(['indLine/ total Lines: ',num2str(indLine), '/', num2str(iLine), ' RF phase: ', num2str(rf_spoiling_phase)]);
        disp([ ' adc phase: ', num2str(adc_phaseOffset{indLine} )]);
    end
  
end

rf_dur = 300;
rf = mr.makeBlockPulse(alpha*pi/180,sys,'Duration',rf_dur*1e-6, 'use', 'excitation');


if add_rfdelay
    rf_delay = mr.makeDelay(rf.deadTime*rfdelay_ratio);
end


%% Set the timing
TR=3.9e-3;
TE = 1.79e-3;
% The Ernst angle
% TE_kernel
if add_rfdelay
    te_kernel = (mr.calcDuration(rf)-rf.delay)/2  + mr.calcDuration(rf_delay) + mr.calcDuration(gzPre) + gz.riseTime+gz.flatTime/2;
else
    te_kernel = (mr.calcDuration(rf)-rf.delay)/2  + mr.calcDuration(gxPre) + mr.calcDuration(gx)/2;
end

% 
if te_kernel<TE
    TE_case = 2;
    switch TE_case
        case 1
            % TE = delay+gxPre+half gx
        delay_te=ceil((TE - mr.calcDuration(gxPre)...
            - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
        case 2
            % TE = RF2+delay+gxPre+half gx
           delay_te=ceil((TE - te_kernel)/seq.gradRasterTime)*seq.gradRasterTime;

        otherwise
            delay_te = abs(rf1.deadTime - mr.calcDuration(rf_delay));
    end
    disp(['TE kernel: ', num2str(te_kernel)])
    disp('TE kernel length < TE ')
    disp(['add delay_te: ', num2str(delay_te)])
    add_delay_te = true;

else
    disp(['TE kernel length is already: ', num2str(te_kernel)])
    disp(['TE', num2str(TE)])
    add_delay_te = false;
    delay_te = 0;
end
% 
% TR_kernel
if add_rfdelay
    tr_kernel = mr.calcDuration(rf) + mr.calcDuration(rf_delay) + delay_te + mr.calcDuration(gz_1)+ mr.calcDuration(gz_2);
else
    tr_kernel = mr.calcDuration(rf) + delay_te + mr.calcDuration(gz_1)+ mr.calcDuration(gz_2);
end

% %
if tr_kernel<TR
    delay_tr=ceil((TR - tr_kernel)/seq.gradRasterTime)*seq.gradRasterTime;
    disp(['TR kernel: ', num2str(tr_kernel)])
    disp('TR kernel length < TR ')
    disp(['add delay_tr: ', num2str(delay_tr)])
    add_delay_tr = true;
else
    delay_tr = 0;
    disp(['TR kernel length is already: ', num2str(tr_kernel)])
    disp(['TR', num2str(TR)])
    add_delay_tr = false;
end


%% 3d radial trajectory


switch traj_design
    case 0
        disp('original traj')
        [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_original(nShot, nSeg, flagSelfNav);
    case 1
        [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_poletopole(nShot, nSeg, flagSelfNav);
    case 2
        [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_continuous(nShot, nSeg, flagSelfNav);
    case 3
        disp('-----------------unishuffle traj-------------------')
        [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_uniform_shuffled(nShot, nSeg, flagSelfNav);
end


%
thetas = [azimuthalAngle',polarAngle'];
disp('thetas size');
disp(size(thetas));
% Loop over phase encodes and define sequence blocks
clear seq;
seq=mr.Sequence(sys);           % Create a new sequence object

%



for indLine = 1:iLine
        polar  = polarAngle(indLine);
        azimut = azimuthalAngle(indLine);
        [gx1,gy1,gz1] = RotateGzWrtZaxis_V02(polar,azimut,gz_1,sys);
        [gx2,gy2,gz2]  = RotateGzWrtZaxis_V02(polar,azimut,gz_2,sys);

        %%%++++++++++++++++++++++++++
        rf.phaseOffset = rf_spoiling_phase_rad{indLine};
        seq.addBlock(rf)
 
        if add_rfdelay
            seq.addBlock(rf_delay);
        end
     
        if add_delay_te
            seq.addBlock(mr.makeDelay(delay_te));
        end
        adc  = mr.makeAdc(Nx,'Duration',gz.flatTime,'Delay',gz.riseTime, ...
            'system',sys, 'phaseOffset', adc_phaseOffset{indLine});
        adc.delay         = adc.delay+mr.calcDuration(gzPre);
        seq.addBlock(gx1,gy1,gz1,adc);
        seq.addBlock(gx2,gy2,gz2);
       
        seq.addBlock(mr.makeDelay(delay_tr));
        if mod(indLine, 500)==0
        str = append('Steady state kspace sampling module kSpace=(',num2str(indLine),'/',num2str(iLine),') generated');
        disp(str);
        end
        TR_kernel = seq.duration();

end

% Time checking
[ok, error_report]=seq.checkTiming;


if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
%% plot
plot_seq = true;
if plot_seq
    seq.plot('timeRange', [0 10]*TR);
end

disp('sequence length:')
disp((seq.duration))


%% Compile the sequence
if do_compile
    % Additionaldefinitions
    
    seq.setDefinition('FOV', fov);
    seq.setDefinition('Name', 'GRE');
    
    
    if isfolder(seqfolder)
        disp('The folder exists: ')
        disp(seqfolder)
    else
        mkdir(seqfolder)
        disp('The folder is created:')
        disp(seqfolder)
    end
    
    seq_name = strcat("", seq_note);
    seqpath = fullfile(seqfolder, strcat(seq_name,'.seq') );
    seq.write(seqpath)       % Write to pulseq file
    disp('The seq is written here:')
    disp(seqpath);
end
%%

% k-space trajectory calculation
kspace_traj = seq.calculateKspacePP;
%% Plotting and checking the trajectory
matpath = fullfile(seqfolder, strcat(seq_name,'.mat') );
save('matpath', 'kspace_traj');
disp('The trajectory is written here:')
disp(matpath);



plot_shot = 1;
plot_samples = false;
make_gif = false;
check_traj = 1;
if check_traj
k_trj = kspace_traj;
k_trj = reshape(k_trj, [3, 480, 22, round(iLine/22)]);
end
%%
if plot_shot
    %plot k-spaces
    % figure;
    % plot(kspace_traj(2,:), kspace_traj(3,:), "o")
    figure ('Color', 'White')
    % for f = 1:size(thetas,1)
    % for iShot = 1:round(iLine/22)
     for iShot = 1:5
        
            kkx = squeeze(k_trj(1,1,:, iShot));
            kky = squeeze(k_trj(2,1,:, iShot));
            kkz = squeeze(k_trj(3,1,:, iShot));
     
            plot3(kkx, kky, kkz,  '-o', 'Markersize', 4, LineWidth=2)
            % exportgraphics(gcf,'sp_interleave_yj.gif','Append',true);
            xlim([-1000,1000])
            ylim([-1000,1000])
            zlim([-1000,1000])
            hold on 
            grid on
            
       
        pause(0.005)
    end
end


if plot_samples
  figure ('Color', 'White')
    % for f = 1:size(thetas,1)
    for iShot = 1
        for iSeg = 1:nSeg
           
            kkx = squeeze(k_trj(1,:,iSeg, iShot));
            kky = squeeze(k_trj(2,:,iSeg, iShot));
            kkz = squeeze(k_trj(3,:,iSeg, iShot));
     
            plot3(kkx, kky, kkz,  '-o', 'Markersize', 4, LineWidth=2)
            % exportgraphics(gcf,'sp_interleave_yj.gif','Append',true);
            xlim([-1000,1000])
            ylim([-1000,1000])
            zlim([-1000,1000])
            hold on 
            grid on
            
        % end
        pause(0.5)
        end
    end
end

if make_gif
    filename = "traj.gif";
    for i = 1:size(kspace_traj,2) % Change loop iterations as needed
      plot3(kspace_traj(1,i), kspace_traj(2,i), kspace_traj(3,i), "o");
      hold on
      xlim([-500,500])
      ylim([-500,500])
      lim([-500,500])
      frame = getframe(gcf);
      im = frame2im(frame);
      [imind, cm] = rgb2ind(im, 256);
      if iff == 1
            imwrite(imind,cm,filename,'gif','LoopCount',inf,'DelayTime',0.2);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
        end
    end
end

clear;

