function [ model ] = choose_cb_model(name,pattern, yoyoConst, netrConst)
    % choose_cb_model
    % Generates theory model parameters for the competitive binding model
    % Puts the parameter values in an n-dimensional matrix, where n is the
    % lenght of a rule
    %
    %     Args:
    %         name: model name
    %     Returns:
    %         model: structure with netropsin and yoyo parameters
    % 
    %     Example:
    %         [ model ] = choose_cb_model()
    %
    %     TODO: add possibility of different models, that were tested in the paper
    %
    
    if nargin < 1
        name = 'literature';
    end

    if nargin < 2
        pattern = 'binding_constant_rules.txt';
    end
    
    if nargin < 4
        yoyoConst = 26;
        netrConst = 0.4;
    end
    
    switch name
        case 'literature'
            % import the binding constant rules based on the model name
            import CBT.Hca.Import.import_binding_constant_rules;
            [bindingConstantNames,bindingConstantVals] = import_binding_constant_rules(name);

            % lengths of small sequences
            seqSpecLen = length(bindingConstantNames(1,:));

            % binding constant matrix, used in CB generation
            bindingConstantsMatSize = repmat(4, [1, seqSpecLen]);
            bindingConstantsMat = NaN(bindingConstantsMatSize);

            % number of binding constant rules
            numRules = size(bindingConstantVals, 1);

            % convert vector to int
            bitsmartTranslationArr = uint8(pow2(seqSpecLen-1:-1:0));

            for ruleNum=1:numRules
                vect_uint8 = bitsmartTranslationArr(nt2int(bindingConstantNames(ruleNum,:)));
                mat_logical = logical(rem(floor(double(vect_uint8(:))*pow2(1 - seqSpecLen:0)),2));
                idxs = mat2cell(mat_logical, ones([1, size(mat_logical, 1)]), 4);
                bindingConstantsMat(idxs{:}) = bindingConstantVals(ruleNum);
            end

            % multiply the constants by parameters from th epaper
            model.netropsinBindingConstant = 0.4*bindingConstantsMat./1E6;
            model.yoyoBindingConstant = 26;
            model.pattern = nt2int('CB');
        case 'enzyme'
            if nargin < 2
                pattern = 'TCGA';
            end
            
            model.pattern = nt2int(pattern);
        case 'custom'
     % import the binding constant rules based on the model name
            import CBT.Hca.Import.import_binding_constant_rules;
            [bindingConstantNames,bindingConstantVals] = import_binding_constant_rules(name,pattern);

            % lengths of small sequences
            seqSpecLen = length(bindingConstantNames{1});

            % binding constant matrix, used in CB generation
            bindingConstantsMatSize = repmat(4, [1, seqSpecLen]);
            bindingConstantsMat = NaN(bindingConstantsMatSize);

            % number of binding constant rules
            numRules = size(bindingConstantVals, 1);

            % convert vector to int
            bitsmartTranslationArr = uint8(pow2(seqSpecLen-1:-1:0));

            for ruleNum=1:numRules
                vect_uint8 = bitsmartTranslationArr(nt2int(bindingConstantNames{ruleNum}));
                mat_logical = logical(rem(floor(double(vect_uint8(:))*pow2(1 - seqSpecLen:1 - seqSpecLen+4-1)),2));
                idxs = mat2cell(mat_logical, ones([1, size(mat_logical, 1)]), 4);
                bindingConstantsMat(idxs{:}) = bindingConstantVals(ruleNum);
            end

            % multiply the constants by parameters from the paper
            model.netropsinBindingConstant = netrConst*bindingConstantsMat./1E6;
            model.yoyoBindingConstant = yoyoConst;
%             model.pattern = nt2int('CB');
        otherwise
            model.pattern = '';
    end
    
end

