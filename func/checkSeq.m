clc;clear;
addpath(genpath('/Users/cag/Documents/forclone/pulseq_v15/pulseq'));
%%
grad_mode = 'Fast';
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

seq=mr.Sequence(sys);

seq.read(strcat('/Users/cag/Documents/forclone/pulseq4mreye/mreye2p0/', ...
    'sub4_gre_pre_oriTraj.seq'));
%% Access block 262
TR = 8.01e-3;
 seq.plot('timeRange', [0 3]*TR);
%%
dur = seq.duration();
%% Display the block's details
disp(block_262);
% Check the RF pulse properties in block 262 (if any)
rf_pulse = block_262.rf;
disp(rf_pulse);

% Check gradient waveforms in block 262 (if any)
gradient_waveform = block_262.gradients;
disp(gradient_waveform);
%%
% k-space trajectory calculation
kspace_traj = seq.calculateKspacePP;
%%
plot_shot = 0;
plot_samples = 1;
make_gif = false;
k_trj = kspace_traj;
nSeg=22;
nLine=419;
N=480;
k_trj = reshape(k_trj, [3, 480, nSeg, nLine]);
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
           for iSample = 1:N
                kkx = squeeze(k_trj(1,iSample,iSeg, iShot));
                kky = squeeze(k_trj(2,iSample,iSeg, iShot));
                kkz = squeeze(k_trj(3,iSample,iSeg, iShot));
         
                plot3(kkx, kky, kkz,  '-o', 'Markersize', 4, LineWidth=2)
                % exportgraphics(gcf,'sp_interleave_yj.gif','Append',true);
                xlim([-1000,1000])
                ylim([-1000,1000])
                zlim([-1000,1000])
                hold on 
                grid on
                
            % end
            pause(0.05)
           end
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

