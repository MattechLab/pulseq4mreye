clc;clear;
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

seq.read(strcat('/Users/cag/Documents/forclone/pulseq4mreye/debug_0402/seqs_v14/', ...
    'sub5_main_libre_rfsp_qua117_adc_ph.seq'));
%% Access block 262
TR = 8.01e-3;
 seq.plot('timeRange', [40000 40010]*TR);

%% Display the block's details
disp(block_262);
% Check the RF pulse properties in block 262 (if any)
rf_pulse = block_262.rf;
disp(rf_pulse);

% Check gradient waveforms in block 262 (if any)
gradient_waveform = block_262.gradients;
disp(gradient_waveform);