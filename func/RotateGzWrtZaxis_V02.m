function [gx,gy,gz] = RotateGzWrtZaxis_V02(polar,azimut,gz_1,sys)
% Author: Nils
% rotate the input gradient around x by polar and z axis by azimuth angle

    % Rotate gz_1 Around the X-Axis
    gz_rot_polar = mr.rotate('x',polar,gz_1);
    
    % % Rotate the Resulting Gradients Around the Z-Axis
   
    array = {};
    count = 1; 
    for k = 1:size(gz_rot_polar,2)
        gz_rot_azi = mr.rotate('z',azimut,gz_rot_polar{k});
        for p = 1:size(gz_rot_azi,2)
            array{count} = gz_rot_azi{p};
            count = count+1;
        end
    end
    lol = array{1};

    prephase_min_dur = lol.shape_dur;
    gx = mr.makeTrapezoid('x','Area',0,"duration", prephase_min_dur,'system',sys);
    gy = mr.makeTrapezoid('y','Area',0,"duration", prephase_min_dur,'system',sys);
    gz = mr.makeTrapezoid('z','Area',0,"duration", prephase_min_dur,'system',sys);

    Nb = size(array,2);     

    for k = 1:Nb
        tmp = array{k};
        str_chan = tmp.channel;
    
        if str_chan=='x'
            gx = tmp;
        end
        if str_chan=='y'
            gy = tmp;
        end
        if str_chan=='z'
            gz = tmp;
        end
    end

end


