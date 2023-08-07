%% function bruker_getAcqPvVersion(Acqp, format)
%   get ParaVision version with which the raw data has been acquired
%
% Input:
%   Acqp:   An acqp struct as generated by the function 
%           readBrukerParamFile('path/acqp')
%   format: (optional) One of the following format specifiers
%           'string': return the complete ParaVision / ParaVision 360
%                     name string
%           'is360' : return boolean whether the raw data was acquired with 
%                     ParaVision or with ParaVision 360.
%           'major' : return major version number (e.g. 5, 6, ... for
%                     ParaVision or 1, 2, ... for ParaVision 360). Always
%                     interpret together with 'is360'.
%
% Output:
%   acqPvVersion: ParaVision version with which the raw data has been 
%                 acquired, according to format specifier.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2019 - 2020
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function acqPvVersion = bruker_getAcqPvVersion(Acqp, format)
    if ~isfield(Acqp,'ACQ_sw_version')
        error('Parameter ACQ_sw_version is missing');
    end
    if ~exist('format', 'var')
        format='string';
    end

    is360 = ~isempty(regexpi(Acqp.ACQ_sw_version, '^PV.360'));
    
    switch format
        case 'is360'
            acqPvVersion = is360;
        case 'major'
            if is360
                [d1, d2, d3, d4, token] = regexpi(Acqp.ACQ_sw_version, ...
                    '^PV.360\.(\d+)\.', 'once');
            else
                [d1, d2, d3, d4, token] = regexpi(Acqp.ACQ_sw_version, ...
                    '^PV.(\d+)\.', 'once');
            end
            acqPvVersion = str2double(token);
        otherwise
            acqPvVersion = Acqp.ACQ_sw_version;
    end
end