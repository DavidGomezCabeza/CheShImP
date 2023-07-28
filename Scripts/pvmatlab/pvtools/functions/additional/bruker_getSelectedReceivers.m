function numSelectedReceivers = bruker_getSelectedReceivers(Acqp, varargin)
% numSelectedReceivers = bruker_getSelectedReceivers(Acqp, varargin)
% Outputs: 
%   numSelectedReceiversNumber (integer): number of selected receivers
%
% Inputs: 
%   Acqp (struct): the Acqp-struct of RawDataObject or readBrukerParamFile
%
% Optional inputs:
%   jobIndex (integer): job index [0..N-1]
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2012
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if length(varargin) == 1
        if isfield(Acqp, 'ACQ_ReceiverSelectPerChan')
            jobIndex=varargin{1};
            jobChan = Acqp.ACQ_jobs(8,jobIndex+1);
            numSelectedReceivers = 0;
            for recNum=1:size(Acqp.ACQ_ReceiverSelectPerChan,2)
                if strcmpi(Acqp.ACQ_ReceiverSelectPerChan(jobChan{1,1}, recNum, :), 'Yes');
                    numSelectedReceivers=numSelectedReceivers+1;
                end
            end
        else
            error('receive selection per job not supported by dataset. please remove job index from function call.')
        end         
    elseif isempty(varargin)
        if strcmp(Acqp.ACQ_experiment_mode,'ParallelExperiment')
            if isfield(Acqp, 'GO_ReceiverSelect')
                if ischar(Acqp.GO_ReceiverSelect)
                    numSelectedReceivers=0;
                    for i=1:size(Acqp.GO_ReceiverSelect,1)
                        if strcmpi(Acqp.GO_ReceiverSelect(i,:), 'Yes');
                            numSelectedReceivers=numSelectedReceivers+1;
                        end
                    end
                else
                    numSelectedReceivers=sum(strcmp(Acqp.GO_ReceiverSelect,'Yes'));
                end
            elseif isfield(Acqp, 'ACQ_ReceiverSelect')
                if ischar(Acqp.ACQ_ReceiverSelect)
                    numSelectedReceivers=0;
                    for i=1:size(Acqp.ACQ_ReceiverSelect,1)
                        if strcmpi(Acqp.ACQ_ReceiverSelect(i,:), 'Yes');
                            numSelectedReceivers=numSelectedReceivers+1;
                        end
                    end
                else
                    numSelectedReceivers=sum(strcmp(Acqp.ACQ_ReceiverSelect,'Yes'));
                end
            else
                disp('No information about active channels available.');
                numSelectedReceivers=1;
                disp('The number of channels is set to 1 ! But at this point the only effect is a bad matrixsize.');
                disp('Later you can change the size yourself.');
            end
        else
            numSelectedReceivers=1;
        end
    else 
        error('too many input arguments!')
    end
    