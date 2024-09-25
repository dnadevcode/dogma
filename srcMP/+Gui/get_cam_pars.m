function [nmbp,nmpx, nmpsf] = get_cam_pars()
    prompt = {'Camera (nm/px)','Point spread function (nm)', 'Extension (nm/bp)'};
    title = 'Camera/experiment params';
    dims = [1 100];
    definput = {'130','300','0.3'};
    answer = inputdlg(prompt,title,dims,definput);

    nmpx  = str2double(answer{1}); % camera nm/px
    nmpsf = str2double(answer{2}); % psf in nanometers
    nmbp  = str2double(answer{3}); % how many delta's should we take

end

