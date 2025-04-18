%% generate k-space trajectory with pulseq:
% please go to the file: /Users/cag/Documents/forclone/pulseq4mreye/debug_0407
% for saving the traj_kspace
% and go to the function for traj for saving the traj_func


traj_kspace = load('/Users/cag/Documents/forclone/pulseq4mreye/func/traj_kspace.mat');
traj_func = load('/Users/cag/Documents/forclone/pulseq4mreye/func/traj_func.mat');
fields = fieldnames(traj_kspace);
traj_kspace = traj_kspace.(fields{1});

fields = fieldnames(traj_func);
traj_func = traj_func.(fields{1});
traj_func = permute(traj_func, [4,1,2,3]);
%%
plot_k_space_shot=1;
if plot_k_space_shot
    %plot k-spaces
    % figure;
    % plot(kspace_traj(2,:), kspace_traj(3,:), "o")
    figure ('Color', 'White')
    % for f = 1:size(thetas,1)
    % for iShot = 1:round(iLine/22)
     for iShot = 1:5
        
            kkx = squeeze(traj_kspace(1,1,:, iShot));
            kky = squeeze(traj_kspace(2,1,:, iShot));
            kkz = squeeze(traj_kspace(3,1,:, iShot));
     
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
%%
plot_traj_func_shot=1;
if plot_traj_func_shot
    %plot k-spaces
    % figure;
    % plot(kspace_traj(2,:), kspace_traj(3,:), "o")
    
    figure ('Color', 'White')
    % for f = 1:size(thetas,1)
    % for iShot = 1:round(iLine/22)
     for iShot = 1:5
        
            kkx = squeeze(traj_func(1,1,:, iShot));
            kky = squeeze(traj_func(2,1,:, iShot));
            kkz = squeeze(traj_func(3,1,:, iShot));
     
            plot3(kkx, kky, kkz,  '-o', 'Markersize', 4, LineWidth=2)
            % exportgraphics(gcf,'sp_interleave_yj.gif','Append',true);
            xlim([-0.5,0.5])
            ylim([-0.5,0.5])
            zlim([-0.5,0.5])
            hold on 
            grid on
            
       
        pause(0.005)
    end
end
%%
R =   max(k_trj(:));
t_monalisa = permute(Traj3D, [4,1,2,3]);

t_pulseq = k_trj/R*0.5;

tol = 1e-4;  % set your preferred tolerance
% Compute the difference
diff = abs(squeeze(t_monalisa(:,:,:)) - squeeze(t_pulseq(:,:,:)));
% Check if all differences are below tolerance
if all(diff(:) < tol)
    disp('✅ All elements are close within the specified tolerance.');
else
    disp('❌ Some elements differ beyond the tolerance.');
    
    % Optionally, show where and how much
    [i, j] = find(diff >= tol);
    for idx = 1:length(i)
        fprintf('❌Mismatch at (%d, %d):  -t = %g, ref = %g, diff = %g\n', ...
            i(idx), j(idx), ...
            t_monalisa(1,i(idx),j(idx)), ...
            t_pulseq(1,i(idx),j(idx)), ...
            diff(i(idx),j(idx)));
    end
end