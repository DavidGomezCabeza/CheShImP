function [output_args, outFlag] = bruker_addFlag(input_args, name, default)
% [output_varargs, outFlag] = bruker_addParamValue(input_varargs, name, default)
% This is a helper function for the input parser to check if a given flag is
% present in the input.
% 
% The input_varargs cell array is searched for the presence of a string given
% as argument 'name' (the name of the flag). If this flag is present, the input
% cell array is cleared of that flag and returned as the output cell array
% 'output_args', and for output 'outFlag', true is returned.
% If the flag is not present, input_args are returned unmodified and for
% output 'outFlag', the value for 'default' is returned.
%
% Example: [output_args, myFlag] = bruker_addParamValue(varargin, 'myFlag', false)
% with varargin={v1, 'MyParam', [0, 1, 2], 'myFlag'}
% output: 
%   output_args={v1, 'MyParam', [0, 1, 2]}
%   outFlag = true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2020
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create mutable copy
    output_args=input_args;

    foundName = strcmpi(output_args, name);

    if 1 == sum(foundName)
        pos = find(1 == foundName);
        outFlag = true;
        % clear flag from cell array
        output_args(pos) = [];
        clear pos;
    else
        if islogical(default)
            % flag not found, so return valid default
            outFlag = default;
        else
            % flag not found, and no valid default given
            outFlag = false;
        end
    end

end
