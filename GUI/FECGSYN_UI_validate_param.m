function [param_out, valid ] = FECGSYN_UI_validate_param( param )
%VALIDATE_PARAM Makes sure param is of the correct format
%   Makes sure all fields in param are valid and adds default values to
%   the fields that are not defined
%
% Mohsan Alvi (mohsan.alvi@eng.ox.ac.uk) - July 2014

try
    disp('TODO: validation of input');
    %% Add default values where missing
    param_out = FECGSYN_UI_add_default_params(param);

    
    
    valid = 1;
    
catch
    param_out = param;
    valid = 0;
    
end


end

