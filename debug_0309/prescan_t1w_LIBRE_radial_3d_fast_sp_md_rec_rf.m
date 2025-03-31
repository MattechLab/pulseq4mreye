%% This file intergrate Nils' phyllotaxis calculation
% The way to calculate timing will be different
% adapted by Yiwei Jia
% LIBRE GRE RADIAL PHYLLOTAXIS

clear;clc;close all;
addpath(genpath('/Users/cag/Documents/forclone/pulseq'))
%
grad_mode = 'Fast';
checking = true;
seq_note = 'prescan_LIBRE_radial_3d_fast_sp_md_rec_rf';
add_rfdelay = true;
rfdelay_ratio = 1;
do_compile = 1;
methodrf = 3; 
%
%%%%%%%%% QQ1: How to determine the param here, 
% should they align with the hardware limit?%%%%%%%%%%%%%%%%
switch grad_mode
    case 'Fast'
        max_grad = 26;      % Max gradient strength [mT/m]
        max_slew = 120;  % Maximum slew rate [mT/m/ms] orginial 120
    case 'Normal'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 120;     % Maximum slew rate [mT/m/ms]
    case 'Whisper'
        max_grad = 12;      % Max gradient strength [mT/m]
        max_slew = 40;      % Maximum slew rate [mT/m/ms]
    case 'Performance'
        max_grad = 40;      % Max gradient strength [mT/m]
        max_slew = 160;  % Maximum slew rate [mT/m/ms]
end



sys = mr.opts('MaxGrad',max_grad,'GradUnit','mT/m',...
    'MaxSlew',max_slew,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6, 'B0', 2.89);

seq=mr.Sequence(sys);           % Create a new sequence object

%% Define the parameters
MatrixSize = 480; %
N = 480; 
% 1/bw: 
bw_px = 496; %Hz/px
adc_dur = 1e6/bw_px;
adc_dwell = round(adc_dur/N/(sys.rfRasterTime*1e6))*(sys.rfRasterTime*1e6);
adc_dur = adc_dwell*N;
%make sure the adc_dwell is valid/integer
% ===Modified end=============================

%This should be the matrix size (matrix size 240, 0.5mm voxel size, 120mm)

% Nx is not necessarily the value in the protocol, should adapt it to the
% duration, in this case: 2016/4 = 504 ??
% But is it legal to just set it as 240?
% Siemens scanner 120mm in protocol but 240mm in raw data?

res = [0.5 0.5 0.5];                        % Define resolution, voxel size (mm)

fov=[MatrixSize*res(1)*1e-3 MatrixSize*res(2)*1e-3 MatrixSize*res(3)*1e-3];     % Define FOV, isotropic, unit[m]
if checking
disp(['fov: ', num2str(fov)])
end

deltak=1 ./ fov;


Nx      = 1e3*fov(1)/(res(2)); 
Ny      = 1e3*fov(2)/(res(2)); 
Nz      = 1e3*fov(3)/res(3);  

% old gradients
% roDuration= ceil((adc_dur*1e-6)/seq.gradRasterTime)*seq.gradRasterTime;         % ADC duration??
% gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak(1),'FlatTime',roDuration,'system',sys);
% gy_rad = mr.makeTrapezoid('y','FlatArea',Nx*deltak(2),'FlatTime',roDuration,'system',sys);
% gz_rad = mr.makeTrapezoid('z','FlatArea',Nx*deltak(3),'FlatTime',roDuration,'system',sys);
% %
% adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
% adc_delay = adc.delay;
% 
% %prephaser: shift to the center of k-space
% gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);
% gyPre = mr.makeTrapezoid('y','Area',-gx.area/2,'system',sys);
% gzPre = mr.makeTrapezoid('z','Area',-gx.area/2,'system',sys);
% 
% % The following two were not called afterwards
% % phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
% % peScales=phaseAreas/gyPre.area;
% 
% % gradient spoiling, the area of spoiler is adjustable...
% gxSpoil=mr.makeTrapezoid('x','Area',0.5*Nx*deltak(1),'system',sys);
% gySpoil=mr.makeTrapezoid('y','Area',0.5*Nx*deltak(2),'system',sys);
% gzSpoil=mr.makeTrapezoid('z','Area',0.5*Nx*deltak(3),'system',sys);

% new gradient from Nils
roDuration= ceil((adc_dur*1e-6)/seq.gradRasterTime)*seq.gradRasterTime;
gz               = mr.makeTrapezoid('z','FlatArea',Nx*deltak(1),'FlatTime',roDuration,'system',sys);
adc              = mr.makeAdc(Nx,'Duration',gz.flatTime,'Delay',gz.riseTime,'system',sys);
gzPre            = mr.makeTrapezoid('z','Area',-gz.area/2,'system',sys);
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


%%RF objects
% Creating the RF object
% How to determine the duration of rf pulse?
alpha = 8;                               % Define flip angle
df = -500;
tao = 1100; %us

phase_accum1 = 1e-6*(0:1:tao-1)*df*2*pi;
b1 = (ones(1,tao));
b1 = b1.*exp(1i*phase_accum1);

phase_accum2 = 1e-6*(0:1:tao-1)*df*2*pi+2*pi*df*tao*1e-6;
b2 = (ones(1,tao));
b2 = b2.*exp(1i*phase_accum2);
if methodrf ==1
    [rf1] = mr.makeArbitraryRf([b1] , alpha*pi/180, sys);
    [rf2] = mr.makeArbitraryRf([b2] , alpha*pi/180, sys);
elseif methodrf == 2 
    rfinterval = complex(zeros(1, 20),0);
    rf = mr.makeArbitraryRf([b1,  rfinterval, b2], alpha*pi/180, sys);
    rf_libre_amp_norm = mr.makeBlockPulse(alpha*pi/180,sys,'Duration',tao*1e-6);
    rf.signal = rf.signal/(max(abs(rf.signal)))*max(abs(rf_libre_amp_norm.signal));
else
    rf_dur = 300;
    rf = mr.makeBlockPulse(alpha*pi/180,sys,'Duration',rf_dur*1e-6);
 
end


if add_rfdelay
    rf_delay = mr.makeDelay(rf.deadTime*rfdelay_ratio);
end
%%
figure;
plot(phase_accum2)
hold on;
plot(phase_accum1);
grid on;


%% Set the timing
TR=8.01e-3;
TE = 3.05e-3;
% 
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
flagSelfNav = 1;
nSeg = 22;
nShot = 419;

nShot_plot = 2;
traj_design = 0;    % Trajectory design: 0 = original, 1 = pole-to-pole, 2 = continuous

switch traj_design
    case 0
        [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_original(nShot, nSeg, flagSelfNav);
    case 1
        [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_poletopole(nShot, nSeg, flagSelfNav);
    case 2
        [polarAngle, azimuthalAngle, vx, vy, vz] = phyllotaxis3D_continuous(nShot, nSeg, flagSelfNav);
end


nLine = size(polarAngle,2);
%
thetas = [azimuthalAngle',polarAngle'];
disp('thetas size');
disp(size(thetas));
% Loop over phase encodes and define sequence blocks
clear seq;
seq=mr.Sequence(sys);           % Create a new sequence object

%

% for indRadial = 1:nLine
nRadial = nLine;
for indRadial = 1:nRadial
        polar  = polarAngle(indRadial);
        azimut = azimuthalAngle(indRadial);
        [gx1,gy1,gz1] = RotateGzWrtZaxis_V02(polar,azimut,gz_1,sys);
        [gx2,gy2,gz2]  = RotateGzWrtZaxis_V02(polar,azimut,gz_2,sys);

        %%%++++++++++++++++++++++++++
        seq.addBlock(rf)
        if add_rfdelay
            seq.addBlock(rf_delay);
        end
     
        if add_delay_te
            seq.addBlock(mr.makeDelay(delay_te));
        end

        seq.addBlock(gx1,gy1,gz1,adc);
        seq.addBlock(gx2,gy2,gz2);
       
        seq.addBlock(mr.makeDelay(delay_tr));

        str = append('Steady state kspace sampling module kSpace=(',num2str(indRadial),'/',num2str(nLine),') generated');
        disp(str);
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
%%plot
plot_seq = true;
if plot_seq
    seq.plot('timeRange', [0 1]*TR);
end

disp('sequence length:')
disp((seq.duration))
%%

% k-space trajectory calculation
kspace_traj = seq.calculateKspacePP;

%% Compile the sequence
if do_compile
    seq.setDefinition('FOV', fov);
    seq.setDefinition('Name', 'GRE');
    
    seqfolder = "/Users/cag/Documents/forclone/pulseq_exercise/seqs_0311";
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

%% Plotting and checking the trajectory

%
plot_shot = false;
plot_samples = true;
make_gif = false;
k_trj = kspace_traj;
k_trj = reshape(k_trj, [3, 480, 22, round(nRadial/22)]);
%%
if plot_shot
    %plot k-spaces
    % figure;
    % plot(kspace_traj(2,:), kspace_traj(3,:), "o")
    figure ('Color', 'White')
    % for f = 1:size(thetas,1)
    for iShot = 1:nRadial
        % for iSeg = 1:nSeg
           
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
            
        % end
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


