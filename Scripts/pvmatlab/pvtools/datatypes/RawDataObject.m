classdef RawDataObject < BrukerDataSuperclass
% RawDataObject
% this class reads and stores raw-MRI-Objects
% you can read ACQP-files, method-files, fid-files and job-files
% 
% [...] means: this is an optional argument
%
%
% Constructor:
% 3 versions to create an object:
% 1.) generate an empty object: rawObj = RawDataObject(); 
% 2.) if you use a standard PraVision pathstructure, use this command to create all data: 
%       rawObj = RawDataObject('path_to_PV_data_directory', ['dataPrecision', yourPrecisionString], ...
%                           ['specified_NRs', your_NR_array], ['specified_Jobs', your_Jobs_array], ... 
%                           ['noFidRead'], ['fid_proc']);
%       -> Requires: standard ParaVision pathstructur
%       -> Does: sets all paths to standard ParaVision-paths and reads the acqp file and the fid-file 
%                (reads fidfile only when 'noFidRead' is NOT set)
%       -> Try: read Method-File
% 3.) use the this command to set all file-paths manually and to create all data:
%       obj=RawDataObject( 'path_to_acqpFile', 'path_to_fidFile', 'path_to_methodFile', ['dataPrecision', yourPrecisionString], ...
%                           ['specified_NRs', your_NR_array], ['specified_Jobs', your_Jobs_array], ... 
%                           ['noFidRead'], ['fid_proc']);
%       -> Does: sets all paths to standard PV-paths and reads the acqp file and the fid-file 
%                (reads fidfile only when 'noFidRead' is NOT set)
%       -> Try: read Method-File
%
%
% Parameters
%
%    Required
%      'path_to_PV_data_directory' - (2.) 
%           Type: string, 
%           The path to your Paravision data directory with the fid and acqp file e.g. '../Testfiles/10'
%      'path_to_acqpFile', 'path_to_fidFile', 'path_to_methodFile' - (3.) 
%           Types: string
%           Pathes incl. dataname: eg. path_to_acqpFile='../Testfiles/10/acqp'
%           If you have only jobs and no fid file, please type the end with /fid, the programm will change it automatically 
%
%    Optional
%      ['dataPrecision', yourPrecisionString] - (2.) & (3.)
%         Type: defined string followed by 'single' or 'double', default is 'double'
%         You can change the precision of the data to single if your memory is not enough for double calculations.
%         e.g. obj = RawDataObject('../Testfiles/10', 'dataPrecision', 'single' );
%      ['specified_NRs', your_NR_array] - (2.) & (3.) 
%          Type: defined string followed by an array, default is [] (= read all NRs)
%          Why use: some MRI-data is very big and it isn't necessary to read the full data everytime,
%          with 'specified_NRs' you can specify some NR-Loops to read.
%          Keep in mind, the first NR loop starts with an 1
%          e.g. obj = RawDataObject('../Testfiles/10', 'specified_NRs', [1 5 7] );
%      ['specified_Jobs', your_Jobs_array] - (2.) & (3.) 
%          Type: defined string followed by an array, default is [] (= read all Jobs)
%          Why use: some MRI-data is very big, and it isn't necessary to read the full data everytime, 
%          with 'specified_Jobs' you can specify some Job-Files to read.
%          Keep in mind, the first Job starts with an 0
%          e.g. obj = RawDataObject('../Testfiles/10', 'specified_Jobs', [0 2] );
%      ['noFidRead'] - (2.) & (3.) 
%          Type: defined string
%          You can use this string to exclude the read-process of the fid-file, you can run it later or only use the parameterfiles.
%          e.g. obj = RawDataObject('../Testfiles/10', 'noFidRead'); will only read the acqp and method file
%      ['fid_proc'] - (2.) & (3.)
%          Type: defined string
%          You can use this string to override reading the native rawdata.job0 and instead read the pre-processed '<EXPNO>/fid', or
%          '<PROCNO>/fid_proc.64' for ParaVision 360 V3.0 and later, respectively. This pre-processed fid, as created by some
%          methods (e.g., CSI, PRESS and STEAM) is created from the raw data by performing phase corrections, averages, combining
%          the receive channels and (for CSI) Fourier transforming in the spatial directions.
%
%
% See also FrameDataObject, CKDataObject, ImageDataObject, readBrukerParamFile,
% readBrukerRaw, convertRawToFrame, convertFrameToCKData,
% readBruker2dseq

%======================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2012 - 2020
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Variables
% for more variables please see BrukerDataClass
properties (SetAccess = public) 
   data={}; %Cellarray, contains (mostly) complex matrices with dimensions (NumberOfChannels, Scansize, NumberOfScansPerChannel), if an fidFile exist it is stored in data{1} 
end


%%  Methods
methods
    %% Constructor
    % The constuctor is descripted in the classdescription   
    function obj = RawDataObject(varargin)
        % Usage and arguments see above

        if isempty(varargin) % -> only create empty object           
            %nothing happens !
        else
            % Init: [Output_args, out ] = bruker_addParamValue( Input_args, name, allowed_string, default)        
            [varargin, dataPrecision] =bruker_addParamValue(varargin, 'dataPrecision', '@(x) strcmp(x, ''single'') || strcmp(x, ''double'')', 'double');
            [varargin, specified_NRs] =bruker_addParamValue(varargin, 'specified_NRs', '@(x) isnumeric(x)', [] );
            [varargin, specified_Jobs] =bruker_addParamValue(varargin, 'specified_Jobs', '@(x) isnumeric(x)', [] );
            [varargin, fid_proc] = bruker_addFlag(varargin, 'fid_proc', false);
            [varargin, noFidRead] = bruker_addFlag(varargin, 'noFidRead', false);
            readfid = ~noFidRead;

            if (length(varargin)==1 || length(varargin)==3) % -> No error!
                %% wirte constructor arguments to variables

                % First step: choose mode 
                if (length(varargin)==1) % -> directory is like ParaVision path structure

                    obj=obj.setDataPath('auto', varargin{1});

                elseif  length(varargin)==3 % -> the directorys are put in manually
                    obj.Filespath.path_to_acqp=varargin{1};
                    obj.Filespath.path_to_fid=varargin{2};
                    obj.Filespath.path_to_method=varargin{3};
                end

                %% start read-processes:
                % read acqpFile:
                [obj.Acqp, obj.HeaderInformation]=readBrukerParamFile(obj.Filespath.path_to_acqp);

                % try reading method-file:
                try
                    [obj.Method, trash]=readBrukerParamFile(obj.Filespath.path_to_method);
                    clear trash;
                catch
                    warning('MATLAB:bruker_note', 'Method-File can''t be read !');
                end 

                if readfid
                    % Read fid-File and jobs
                    args={obj.Acqp, obj.Method, obj.Filespath.path_to_fid, ...
                        'precision', dataPrecision, 'specified_NRs', specified_NRs, ...
                        'specified_Jobs', specified_Jobs};
                    if fid_proc
                        args = adaptArgsFidProc(obj, args);
                    end
                    [ obj.data, addInfo] = readBrukerRaw(args{:});
                    obj.Acqp = bruker_combineStructs( obj.Acqp, addInfo );
                end
            else
                error('Arguments are not correct');
            end

        end
    end

    %% Methods:
        % the function uses parameters of Acqp to read the fid-File
        % use with: e.g. objectname: obj1 -> obj1=obj1.readRaw();
    function obj=readRaw(obj, varargin)
        % reads the fid and job-files: requires paths in Filespath and
        % Acqp, optional: ['dataPrecision', yourPrecisionString], ['specified_NRs', your_NR_array], 
        %                 ['specified_Jobs', your_Jobs_array], ['fid_proc']

        % Init: % [vararg_out, out_value] =bruker_readArguments(vararg_in, paramName, allowed, default)         
        [trash, dataPrecision] =bruker_addParamValue(varargin, 'dataPrecision', '@(x) strcmp(x, ''single'') || strcmp(x, ''double'')', 'double');
        [trash, specified_NRs] =bruker_addParamValue(varargin, 'specified_NRs', '@(x) isnumeric(x)', [] );
        [trash, specified_Jobs] =bruker_addParamValue(varargin, 'specified_Jobs', '@(x) isnumeric(x)', [] );
        [trash, fid_proc] = bruker_addFlag(varargin, 'fid_proc', false);
        clear trash;

        if isempty(obj.Acqp) || isempty(obj.Filespath.path_to_fid)
            error('Acqp or Filespath.path_to_fid are empty !');
        end
        if ~isempty(specified_Jobs)
            obj.Acqp.specified_Jobs=specified_Jobs;
        end
        if ~isempty(specified_NRs)
            obj.Acqp.specified_NRs=specified_NRs;
        end

        args={obj.Acqp, obj.Method, obj.Filespath.path_to_fid, 'precision', dataPrecision, ...
            'specified_NRs', specified_NRs, 'specified_Jobs', specified_Jobs};
        if fid_proc
            args = adaptArgsFidProc(obj, args);
        end
        [ obj.data ] = readBrukerRaw(args{:});
    end

    function obj=load(obj, load_struct)
        % loading an Object from a struct, saved with struct(obj);
        obj=obj.load_superclass(load_struct);
        obj.data=load_struct.data;
    end
end

methods (Access = private)
    function args = adaptArgsFidProc(obj, args)
        args=[args, {'fid_proc'}];
        if bruker_getAcqPvVersion(args{1}, 'is360') && ...
            bruker_getAcqPvVersion(args{1}, 'major') >= 3
            % ParaVision 360 V3.0 and later store the pre-pocessed
            % fid in the procno folder
            pathToFidProc = obj.Filespath.path_to_reco;
            if length(pathToFidProc) < 4 || ~strcmp(pathToFidProc(end-3:end),'reco')
                error(['Please specify the PROCNO folder containing the ' ...
                    'pre-processed fid (no longer stored in the EXPNO folder ' ...
                    'starting from ParaVision 360 V3.0).']);
            end
            pathToFidProc = [pathToFidProc(1:end-4) 'fid_proc.64'];
            args{3} = pathToFidProc;
        end
    end
end

end
