%%
clear;clc;close all;
addpath(genpath('/Users/cag/Documents/forclone/pulseq'))
%
grad_mode = 'Fast';
checking = true;
seq_note = 'fast_slew120_2p2klines_gapRfGx2deadTime';
add_rfdelay = true;
rfdelay_ratio = 2;
do_compile = false;
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
adc_dur = 480*4;
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

%%%% gradients
roDuration= ceil((adc_dur*1e-6)/seq.gradRasterTime)*seq.gradRasterTime;         % ADC duration??
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak(1),'FlatTime',roDuration,'system',sys);
gy_rad = mr.makeTrapezoid('y','FlatArea',Nx*deltak(2),'FlatTime',roDuration,'system',sys);
gz_rad = mr.makeTrapezoid('z','FlatArea',Nx*deltak(3),'FlatTime',roDuration,'system',sys);
%
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
adc_delay = adc.delay;

%prephaser: shift to the center of k-space
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);
gyPre = mr.makeTrapezoid('y','Area',-gx.area/2,'system',sys);
gzPre = mr.makeTrapezoid('z','Area',-gx.area/2,'system',sys);

% The following two were not called afterwards
% phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
% peScales=phaseAreas/gyPre.area;

% gradient spoiling, the area of spoiler is adjustable...
gxSpoil=mr.makeTrapezoid('x','Area',0.5*Nx*deltak(1),'system',sys);
gySpoil=mr.makeTrapezoid('y','Area',0.5*Nx*deltak(2),'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',0.5*Nx*deltak(3),'system',sys);

% RF objects
% Creating the RF object
% How to determine the duration of rf pulse?
alpha = 8;                               % Define flip angle

df = -500;
phase_accum1 = 1e-6*(0:1:1099)*df*2*pi;
b1 = (ones(1,1100));
b1 = b1.*exp(1i*phase_accum1);

phase_accum2 = 1e-6*(0:1:1099)*df*2*pi+2*pi*df*1099*1e-6;
b2 = (ones(1,1100));
b2 = b2.*exp(1i*phase_accum2);

[rf1] = mr.makeArbitraryRf([b1] , alpha*pi/180, sys);
[rf2] = mr.makeArbitraryRf([b2] , alpha*pi/180, sys);

rf_libre_amp_norm = mr.makeBlockPulse(alpha*pi/180,sys,'Duration',1100*1e-6);
rf1.signal = rf1.signal/(max(abs(rf1.signal)))*max(abs(rf_libre_amp_norm.signal));
rf2.signal = rf2.signal/(max(abs(rf2.signal)))*max(abs(rf_libre_amp_norm.signal));
%%%%%%%%%%%%%%%%% QQ6: why do we choose 2e-5?%%%%%%%%%%%%%
rf1_delay = mr.makeDelay(rf1.deadTime);
if add_rfdelay
    rf2_delay = mr.makeDelay(rf2.deadTime*rfdelay_ratio);
end

%% Set the timing
TR=8.01e-3;
TE = 3.62e-3;

% TE_kernel
if add_rfdelay
    te_kernel = mr.calcDuration(rf2)  + mr.calcDuration(rf2_delay) + mr.calcDuration(gxPre) + mr.calcDuration(gx)/2;
else
    te_kernel = mr.calcDuration(rf2)  + mr.calcDuration(gxPre) + mr.calcDuration(gx)/2;
end

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
    disp('TE kernel length < TE ')
    disp(['add delay_te: ', num2str(delay_te)])
    add_delay_te = true;
    
else
    disp(['TE kernel length is already: ', num2str(te_kernel)])
    disp(['TE', num2str(TE)])
    add_delay_te = false;
    delay_te = 0;
end

% TR_kernel
if add_rfdelay
    tr_kernel = mr.calcDuration(rf1) + mr.calcDuration(rf2) + mr.calcDuration(rf2_delay) + delay_te + mr.calcDuration(gxPre)+ mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
else
    tr_kernel = mr.calcDuration(rf1) + mr.calcDuration(rf2) + delay_te + mr.calcDuration(gxPre)+ mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
end

%
if tr_kernel<TR
    delay_tr=ceil((TR - tr_kernel)/seq.gradRasterTime)*seq.gradRasterTime;
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
nShot = 2055;
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
%%
thetas = [azimuthalAngle',polarAngle'];
disp('thetas size');
disp(size(thetas));
% Loop over phase encodes and define sequence blocks
clear seq;
seq=mr.Sequence(sys);           % Create a new sequence object

%%
% for i=1:size(thetas, 1)
for i = 1:220
    for c=1:length(TE)
        % in our case length(TE)=1, since it is a single echo sequence
        % is the scaler correct or not?
        scalers = GradScalers(thetas(i,1), thetas(i,2));
        %%%++++++++++++++++++++++++++
        seq.addBlock(rf1)
        % seq.addBlock(rf1_delay);
        seq.addBlock(rf2)
        if add_rfdelay
            seq.addBlock(rf2_delay);
        end
     
        gxPre.delay = delay_te;
        gyPre.delay = delay_te;
        gzPre.delay = delay_te;

        %%%++++++++++++++++++++++++++
        seq.addBlock(mr.scaleGrad(gxPre, scalers(1)),...
            mr.scaleGrad(gyPre,scalers(2)),...
            mr.scaleGrad(gzPre,scalers(3)));
        
        seq.addBlock(mr.scaleGrad(gx,scalers(1))...
            ,mr.scaleGrad(gy_rad, scalers(2)),...
           mr.scaleGrad(gz_rad, scalers(3)), adc);
       

        seq.addBlock(gxSpoil,gySpoil,gzSpoil);
        % seq.addBlock(mr.makeDelay(delayTR(c)), gxSpoil, gySpoil, gzSpoil);
        seq.addBlock(mr.makeDelay(delay_tr));

        str = append('Steady state kspace sampling module kSpace=(',num2str(i),'/',num2str(nLine),') generated');
        disp(str);
    end
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

plot_seq = true;
if plot_seq
    seq.plot('timeRange', [0 1]*TR);
end
TR_kernel = (seq.duration);
disp('sequence length:')
disp((seq.duration))
%%

% k-space trajectory calculation
kspace_traj = seq.calculateKspacePP;
% Plot the sequence

% figure
figure
plot(kspace_traj(2,:), kspace_traj(3,:), "o")
figure
plot(kspace_traj(1,:), kspace_traj(2,:), "o")

%% Compile the sequence
if do_compile
    seq.setDefinition('FOV', fov);
    seq.setDefinition('Name', 'GRE');
    
    seqfolder = "/Users/cag/Documents/forclone/pulseq_exercise/seqs";
    if isfolder(seqfolder)
        disp('The folder exists: ')
        disp(seqfolder)
    else
        mkdir(seqfolder)
        disp('The folder is created:')
        disp(seqfolder)
    end
    
    seq_name = strcat("YJ_TR8p01_0p5x0p5x0p5_LIBRE_traj_original_", seq_note);
    seqpath = fullfile(seqfolder, strcat(seq_name,'.seq') );
    seq.write(seqpath)       % Write to pulseq file
    disp('The seq is written here:')
    disp(seqpath);
end

%% Plotting and checking the trajectory

%
plot_shot = true;
plot_samples = false;
make_gif = false;

if plot_shot
    %plot k-spaces
    % figure;
    % plot(kspace_traj(2,:), kspace_traj(3,:), "o")
    figure ('Color', 'White')
    % for f = 1:size(thetas,1)
    for f = 1:22
        kinit = (f-1)*nSeg+1;
        kend = (f)*nSeg;
        plot3(kspace_traj(1,kinit), kspace_traj(2,kinit), kspace_traj(3,kinit), "o")
        % exportgraphics(gcf,'sp_interleave_yj.gif','Append',true);
        xlim([-500,500])
        ylim([-500,500])
        zlim([-500,500])
    hold on 
    grid on
    % pause(5e-10)
    end
end


if plot_samples
    grid on
    axis image
    for i = 1:size(kspace_traj,2)
     % for i = 1:200*480
        plot3(kspace_traj(1,i), kspace_traj(2,i), kspace_traj(3,i), "o");
        % exportgraphics(gcf,'testAnimated.gif','Append',true);
        xlim([-500,500])
        ylim([-500,500])
        zlim([-500,500])
        hold on 
        % pause(0.0000000001)
        grid on
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


%% Function zone
function [m_adPolarAngle, m_adAzimuthalAngle, x, y, z] = phyllotaxis3D_continuous(m_lNumberOfFrames, m_lProjectionsPerFrame, flagSelf)
% adapted by Eva Peper, evaspeper@gmail.com, 11/06/2024
    
    NProj = m_lNumberOfFrames * m_lProjectionsPerFrame ; % = shot x segments
    lTotalNumberOfProjections = NProj;
    
    m_adAzimuthalAngle=zeros(1,NProj);
    m_adPolarAngle=zeros(1,NProj);
    
    x = zeros (1, NProj);
    y = zeros (1, NProj);
    z = zeros (1, NProj);
    
    disp(['NProj', num2str(NProj)]);
    disp(['m_lNumberOfFrames', num2str(m_lNumberOfFrames)]);

    if flagSelf
        N = lTotalNumberOfProjections - m_lNumberOfFrames;
    else
        N = lTotalNumberOfProjections ;
    end
    kost = pi/(2*sqrt(N/2)); % continuous phyllotaxis ESP
    
    Gn = (1 + sqrt(5))/2;
    Gn_ang = 2*pi - (2*pi / Gn);
    count = 1;

    for lk = 1:m_lProjectionsPerFrame
        for lFrame = 1:m_lNumberOfFrames
    
            linter = lk + (lFrame-1) * m_lProjectionsPerFrame;
    
            if flagSelf && lk == 1
    
                m_adPolarAngle(linter) = 0;
                m_adAzimuthalAngle(linter) = 0;
    
            else
    
                % continuous phyllotaxis ESP
                if count<=(N/4)
                    m_adPolarAngle(linter) = kost * sqrt(2*count);
                elseif count>(N/4) && count<=(N/2)
                    m_adPolarAngle(linter) =  pi - kost * sqrt(2*(N/2-count)) ;
                elseif count>(N/2) && count<=(N/4*3)
                    m_adPolarAngle(linter) = pi +  kost * sqrt(2*(count-N/2)) ;
                else
                    m_adPolarAngle(linter) = 2*pi - kost * sqrt(2*(N-count)) ;
                end
    
                m_adAzimuthalAngle(linter) = mod ( (count)*Gn_ang, (2*pi) );
                count = count + 1;
    
            end
    
            x(linter)= sin(m_adPolarAngle(linter))*cos(m_adAzimuthalAngle(linter));
            y(linter)= sin(m_adPolarAngle(linter))*sin(m_adAzimuthalAngle(linter));
            z(linter)= cos(m_adPolarAngle(linter));
    
        end
    end
end

%% I have doubt in this function
function res = GradScalers(azimuth, polar)

resz = cos(polar);
resy = sin(polar)*sin(azimuth);
resx = sin(polar)*cos(azimuth);

res = [resx, resy, resz];

end
%%

function [m_adPolarAngle, m_adAzimuthalAngle, x, y, z] = phyllotaxis3D_original(m_lNumberOfFrames, m_lProjectionsPerFrame, flagSelf)
% original phyllotaxis 

NProj = m_lNumberOfFrames * m_lProjectionsPerFrame; % = shot x segments
lTotalNumberOfProjections = NProj;

m_adAzimuthalAngle=zeros(1,NProj);
m_adPolarAngle=zeros(1,NProj);

x = zeros (1, NProj);
y = zeros (1, NProj);
z = zeros (1, NProj);

if flagSelf
    N = lTotalNumberOfProjections - m_lNumberOfFrames; 
else
    N = lTotalNumberOfProjections ; 
end
kost = pi/(2*sqrt(N));

Gn = (1 + sqrt(5))/2;
Gn_ang = 2*pi - (2*pi / Gn);
count = 1;

for lk = 1:m_lProjectionsPerFrame
    for lFrame = 1:m_lNumberOfFrames

        linter = lk + (lFrame-1) * m_lProjectionsPerFrame;

        if flagSelf && lk == 1

            m_adPolarAngle(linter) = 0;
            m_adAzimuthalAngle(linter) = 0;

        else

            % original phyllotaxis 
            m_adPolarAngle(linter) = kost * sqrt(count);
            m_adAzimuthalAngle(linter) = mod ( (count)*Gn_ang, (2*pi) );
            count = count + 1;

        end

        x(linter)= sin(m_adPolarAngle(linter))*cos(m_adAzimuthalAngle(linter));
        y(linter)= sin(m_adPolarAngle(linter))*sin(m_adAzimuthalAngle(linter));
        z(linter)= cos(m_adPolarAngle(linter));

    end
end

end