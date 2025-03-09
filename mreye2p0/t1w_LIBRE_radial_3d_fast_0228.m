%%
clear;clc;close all;
addpath(genpath('/Users/cag/Documents/forclone/pulseq'))
%%
grad_mode = 'Fast';
checking = true;
%
%%%%%%%%% QQ1: How to determine the param here, 
% should they align with the hardware limit?%%%%%%%%%%%%%%%%
switch grad_mode
    case 'Fast'
        max_grad = 26;      % Max gradient strength [mT/m]
        max_slew = 120;  % Maximum slew rate [mT/m/ms]
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
if false
% ( BW = (1/adc.dwell)/Nx, or in our calc. 1e6/adc_dur )
% ===Modified based on LIBRE BW: 496 Hz/Px===
% dur: int num, dwell is multiple of some value
adc_dur = 2016;
Nx = 502; 
end


Nx = 480; 
adc_dur = 480*4;
% ===Modified end=============================

%This should be the matrix size (matrix size 240, 0.5mm voxel size, 120mm)

% Nx is not necessarily the value in the protocol, should adapt it to the
% duration, in this case: 2016/4 = 504 ??
% But is it legal to just set it as 240?
% Siemens scanner 120mm in protocol but 240mm in raw data?

res = [0.5 0.5 0.5];                        % Define resolution

fov=[Nx*res(1)*1e-3 Nx*res(2)*1e-3 Nx*res(3)*1e-3];     % Define FOV, isotropic, unit[m]
if checking
disp(['fov: ', num2str(fov)])
end

deltak=1 ./ fov;



%%%% gradients
roDuration= ceil((adc_dur*1e-6)/seq.gradRasterTime)*seq.gradRasterTime;         % ADC duration??
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak(1),'FlatTime',roDuration,'system',sys);
gy_rad = mr.makeTrapezoid('y','FlatArea',Nx*deltak(2),'FlatTime',roDuration,'system',sys);
gz_rad = mr.makeTrapezoid('z','FlatArea',Nx*deltak(3),'FlatTime',roDuration,'system',sys);

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

df = -300;
phase_accum1 = 1e-6*(0:1:1309)*df*2*pi;
b1 = (ones(1,1310));
b1 = b1.*exp(1i*phase_accum1);

phase_accum2 = 1e-6*(0:1:1309)*df*2*pi+2*pi*df*1310*1e-6;
b2 = (ones(1,1310));
b2 = b2.*exp(1i*phase_accum2);

[rf1] = mr.makeArbitraryRf([b1] , alpha*pi/180, sys);
[rf2] = mr.makeArbitraryRf([b2] , alpha*pi/180, sys);

rf_libre_amp_norm = mr.makeBlockPulse(alpha*pi/180,sys,'Duration',1300*1e-6);
rf1.signal = rf1.signal/(max(abs(rf1.signal)))*max(abs(rf_libre_amp_norm.signal));
rf2.signal = rf2.signal/(max(abs(rf2.signal)))*max(abs(rf_libre_amp_norm.signal));
%%%%%%%%%%%%%%%%% QQ6: why do we choose 2e-5?%%%%%%%%%%%%%
rf_delay = mr.makeDelay(2e-5);


%% Set the timing
TR=8.01e-3;
TE = 3.62e-3;

% for 3d, no gz is counted.
% delayTE=ceil((TE - mr.calcDuration(gxPre)...
%     - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
% delayTR=ceil((TR - mr.calcDuration(gxPre) ...
%     - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
% assert(all(delayTE>=0));
% assert(all(delayTR>=mr.calcDuration(gxSpoil)));

delay_te = abs(rf1.deadTime - mr.calcDuration(rf_delay));
TR_kernel = 7.6e-3;
delay_tr = ( (((TR-TR_kernel)/2)));

delay_te = ceil(delay_te/seq.gradRasterTime)*seq.gradRasterTime;
delay_tr = ceil(delay_tr/seq.gradRasterTime)*seq.gradRasterTime;

%% 3d radial trajectory

nSeg = 22;
nShot = 2055;
nShot_plot = 2;
[m_adPolarAngle, m_adAzimuthalAngle, x, y, z] = phyllotaxis3D_original(nSeg, nShot, 1);
% m_adPolarAngle = m_adPolarAngle(1:nSeg*nShot_plot);
% m_adAzimuthalAngle = m_adAzimuthalAngle(1:nSeg*nShot_plot);
m_adPolarAngle = m_adPolarAngle(1:nSeg*nShot);
m_adAzimuthalAngle = m_adAzimuthalAngle(1:nSeg*nShot);

thetas = [m_adAzimuthalAngle',m_adPolarAngle'];
disp('thetas size');
disp(size(thetas));
% Loop over phase encodes and define sequence blocks
clear seq;
seq=mr.Sequence(sys);           % Create a new sequence object

for i=1:size(thetas, 1)
% for i = 1:2
    for c=1:length(TE)
        % in our case length(TE)=1, since it is a single echo sequence
        scalers = GradScalers(thetas(i,1), thetas(i,2));
        %%%++++++++++++++++++++++++++
        seq.addBlock(rf1)
        %%%++++++++++++++++++++++++++
        seq.addBlock(rf2, rf_delay)

        
        delay_te = abs(rf1.deadTime - mr.calcDuration(rf_delay));
        gxPre.delay = delay_te+delay_tr;
        gyPre.delay = delay_te+delay_tr;
        gzPre.delay = delay_te+delay_tr;

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
    seq.plot('timeRange', [0 2]*TR);
end
TR_kernel = (seq.duration);
disp('sequence length:')
disp((seq.duration))
%%

% k-space trajectory calculation
kspace_traj = seq.calculateKspacePP;
%% Plot the sequence

% figure
figure
plot(kspace_traj(2,:), kspace_traj(3,:), "o")
figure
plot(kspace_traj(1,:), kspace_traj(2,:), "o")

%% Compile the sequence

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

seq_name = "YJ_TR8p01_0p5x0p5x0p5_noGrappa_wocalib_LIBRE_traj_original";
seqpath = fullfile(seqfolder, strcat(seq_name,'.seq') );
seq.write(seqpath)       % Write to pulseq file
disp('The seq is written here:')
disp(seqpath);


%% Plotting and checking the trajectory

%
plot_interleave = true;
plot_lines = false;
make_gif = false;

if plot_interleave
    %plot k-spaces
    % figure;
    % plot(kspace_traj(2,:), kspace_traj(3,:), "o")
    figure ('Color', 'White')
    for f = 1:size(thetas,1)
        kinit = (f-1)*nSeg+1;
        kend = (f)*nSeg;
        plot3(kspace_traj(1,kinit), kspace_traj(2,kinit), kspace_traj(3,kinit), "o")
        exportgraphics(gcf,'sp_interleave.gif','Append',true);
        xlim([-500,500])
        ylim([-500,500])
        zlim([-500,500])
    hold on 
    pause(0.000000005)
    end
end


if plot_lines
    grid on
    axis image
    for i = 1:size(kspace_traj,2)
        plot3(kspace_traj(1,i), kspace_traj(2,i), kspace_traj(3,i), "o");
        exportgraphics(gcf,'testAnimated.gif','Append',true);
        xlim([-500,500])
        ylim([-500,500])
        zlim([-500,500])
        hold on 
        pause(0.0000000001)
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


function res = GradScalers(azimuth, polar)

resz = cos(polar);
resy = sin(polar)*sin(azimuth);
resx = sin(polar)*cos(azimuth);

res = [resx, resy, resz];

end


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