function [bindingConstantNames,bindingConstantVals] = import_binding_constant_rules(name, filename)
    % import_binding_constant_rules
    
    % 
     % A => 1 , C => 2, G => 3,  T(U) => 4

     if nargin < 2
         filename = 'binding_constant_rules.txt';
     end
    
    switch name
        case 'literature'
            % A => 1 , C => 2, G => 3,  T(U) => 4
            fid = fopen(filename); 
            value = textscan(fid,'%4c %f','delimiter','\n'); fclose(fid);     
            
            bindingConstantNames = value{1};
            
            bindingConstantVals = value{2};
        case 'custom'
            % A => 1 , C => 2, G => 3,  T(U) => 4
            values = importdata(filename);
%             fid = fopen(filename); 
%             value = textscan(fid,'%4c %f','delimiter','\n'); fclose(fid);     
            
            bindingConstantNames = values.textdata;
            
            bindingConstantVals = values.data;
              
        otherwise
            error('binding constant rules undefined');
    end

end

