function [ correctChromosome, correctPlace,allowedError, thresh] = choose_correct_place(defaultChromosome, defaultPlace, defaultError, defaultThresh,titleText)
    
    correctChromosome = defaultChromosome;
    correctPlace = defaultPlace;
    allowedError = defaultError;
    thresh  = defaultThresh;
    options.Resize = 'on';
    % sets number of time-frames
    answer = inputdlg(...
        { ...
            'Choose the number (as appears in the plot) of correct chromosome', 'Choose correct place (pixel resolution) along the chromosome', 'Choose how many pixel error to allow','Choose p-value threshold'.....
        }, ... % prompt
        titleText, ... % dialog title
        1, ... % number of lines
        { ...
            num2str(correctChromosome), ...
            num2str(correctPlace), ...
            num2str(defaultError), ...
            num2str(thresh), ...
        },options);
    if ~isempty(answer)
        correctChromosome =  str2double(answer{1});
        correctPlace = str2double(answer{2});
        allowedError = str2double(answer{3});
        thresh = str2double(answer{4});
    end

    isAcceptable = struct;
    isAcceptable.correctChromosome = (correctChromosome >= 0);
    isAcceptable.correctPlace = (correctPlace >= 0);
    isAcceptable.allowedError = (allowedError >= 0);
    isAcceptable.thresh = (thresh >= 0);

    if not(isAcceptable.correctPlace) ||not(isAcceptable.allowedError) 
        correctChromosome = defaultChromosome;
        correctPlace = defaultPlace;
        allowedError = defaultError;
        thresh = defaultThresh;
        warning('Bad input ! Try again!');
    end
end