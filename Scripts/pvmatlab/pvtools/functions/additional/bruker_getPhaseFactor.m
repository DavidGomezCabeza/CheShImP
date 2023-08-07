%% function phaseFactor = bruker_getPhaseFactor(Acqp, [Method])
%   For PV 6 and 7, extract the phase factor from the Acqp structure. For 
%   PV 360 job acquisition, the phase factor is derived in a method 
%   specific way from the PVM parameters.
% 
%   When objects (e.g., slices) are acquired in an interleaved way, the 
%   phase factor determines the number of subsequent phase encoding steps 
%   acquired in a single object, before raw data from the next object is 
%   acquired.
%
% Input:
%   Acqp: An acqp struct as generated by the function
%         readBrukerParamFile('path/acqp')
%
% Optional Input (required for PV 360):
%   Method: A method struct as generated by the function
%           readBrukerParamFile('path/method')
%
% Output:
%   phaseFactor for the acquired raw data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2021
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phaseFactor = bruker_getPhaseFactor(Acqp, Method)

    % ParaVision
    if ~bruker_getAcqPvVersion(Acqp, 'is360')
        if ~isfield(Acqp, 'ACQ_phase_factor')
            error('ACQ_phase_factor required.');
        end
        phaseFactor = Acqp.ACQ_phase_factor;
        return;
    end
    
    % ParaVision 360
    if ~exist('Method','var')
        error('Method parameters are required to create PV 360 frame data.');
    end
    
    if ~bruker_requires({Method}, {'Method'})
        error('Parameter "Method" is missing.');
    end
    
    % For RARE and variants RAREst, RAREVTR, FAIR_RARE, use PVM_RareFactor
    % FISP and B1Map also use PVM_RareFactor
    if ~isempty(strfind(Method.Method, 'RARE')) ...
            || ~isempty(strfind(Method.Method, 'FISP')) ...
            || ~isempty(strfind(Method.Method, 'B1Map'))
        if ~isfield(Method, 'PVM_RareFactor')
            error('Parameter "PVM_RareFactor" is missing.');
        end
        phaseFactor = Method.PVM_RareFactor;
        return;
    end
    
    % For SegFLASH and MDEFT, use SegmentSize
    if ~isempty(strfind(Method.Method, 'SegFLASH')) ...
            || ~isempty(strfind(Method.Method, 'MDEFT'))
        if ~isfield(Method, 'SegmentSize')
            error('Parameter "SegmentSize" is missing.');
        end
        phaseFactor = Method.SegmentSize;
        return;
    end

    % FLASH and its remaining unsegmented variants have an AngioMode, in 
    % which each all 2D phase encodes are acquired in one go
    if ~isempty(strfind(Method.Method, 'FLASH'))
        if ~isfield(Method, 'AngioMode')
            error('Parameter "AngioMode" is missing.');
        end
        if strcmp(Method.AngioMode,'Yes')
            phaseFactor = Method.PVM_EncMatrix(2);
        else
            phaseFactor = 1;
        end
        return;
    end

    % DESS 
    if ~isempty(strfind(Method.Method, 'DESS'))
        if ~isfield(Acqp, 'ACQ_dim')
            error('Parameter "ACQ_dim" is missing.');
        end
        if 2 == Acqp.ACQ_dim
            phaseFactor = Method.PVM_EncMatrix(2);
        else
            phaseFactor = 1;
        end
        return;
    end

    % All other cases
    phaseFactor = 1;
end