function [spectra, abscissaValuesHz, abscissaValuesPpm] = bruker_loadSpectra(procno, basicPhaseCorrection)
%bruker_loadSpectra Load spectra from pre-processed time-domain data.
%
% [spectra, abscissaValuesHz, abscissaValuesPpm] = ...
%     bruker_loadSpectra(procno, [basicPhaseCorrection]);
%
% This function loads preprocessed time-domain data (fid_proc files), 
% generated by the reconstructions of some Bruker spectroscopic methods and 
% spectroscopic imaging methods. The preprocessed time-domain data is then
% Fourier transformed, resulting in a complex spectrum for each voxel in
% the dataset.
%
% Input:
%   procno                procno path of a reconstructed spectroscopic 
%                         dataset or spectroscopic imaging dataset
%
% Optional input:
%   basicPhaseCorrection  if true, a very basic linear phase correction is
%                         performed. Default is false.
%
% Output:
%   spectra               array of complex spectra (nPoints spectral points)
%                         Dimension: (nPoints, dim1, dim2, dim3, nRepetitions),
%                         with spatial dimensions dim1, dim2, and dim3.
%                         For pure spectroscopic methods, the spatial 
%                         dimensions are all equal to 1.
%
%   abscissaValuesHz      the corresponding abscissa value for each spectral 
%                         point, in Hz. The abscissa values in the array 
%                         are in ascending order.
%
%   abscissaValuesPpm     the corresponding abscissa value for each spectral
%                         point, in ppm. The abscissa values in the array 
%                         are in ascending order.
%
% Note:
%   To display the spectra in spectroscopic convention (descending
%   abscissa values), use the 'xdir' property of the axis. For example:
%     
%     plot(abscissaValuesPpm, abs(spectra(:,1,1,1,1)));
%     set(gca, 'xdir','reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2021
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default for optional parameter, if no value is given.
if ~exist('basicPhaseCorrection','var')
    basicPhaseCorrection = false;
end

% Import required parameters.
Reco = readBrukerParamFile([procno filesep 'reco']);

% Import pre-processed time-domain data.
rawObj = RawDataObject(procno, 'fid_proc');

isEPSI = ~isempty(strfind(rawObj.Acqp.ACQ_method,'EPSI'));

nRepetitions = Reco.RecoNumRepetitions;
spatialDims = length(Reco.RECO_size)-1;
spatialSizes = ones(1,3);

if isEPSI
    % In EPSI, the spectral dimension is in between the first two spatial
    % dimensions
    spatialSizes(1) = Reco.RECO_size(1);
    nPoints = Reco.RECO_size(2);
    spatialSizes(2) = Reco.RECO_size(3);
    if 3 == spatialDims
        spatialSizes(3) = Reco.RECO_size(4);
    end
else
    nPoints = Reco.RECO_size(1);
    for dim = 1:spatialDims
        spatialSizes(dim) = Reco.RECO_size(dim+1);
    end
end

if 2 == spatialDims
    spatialSizes(3) = Reco.RecoObjectsPerRepetition;
end
nVoxels = prod(spatialSizes)*nRepetitions;

timeDomain = rawObj.data{1};
if numel(timeDomain) ~= nPoints * nVoxels
    error('Error: reco and fid_proc files are inconsistent.');
end

if isEPSI
    % sort spectral dimension to the front
    timeDomain = reshape(timeDomain,[spatialSizes(1) nPoints spatialSizes(2) spatialSizes(3) nRepetitions]);
    timeDomain = permute(timeDomain,[2,1,3,4,5]);
    timeDomain = reshape(timeDomain,[1 nPoints nVoxels]);
end

% Perform Fourier-transformation on first direcion, with 0 Hz in the
% center.
ftshift=exp(1i*pi*(0:(nPoints-1))); % 1 -1 1 -1 ...
spectra = timeDomain .* repmat(ftshift, [1, 1, nVoxels]);
spectra = fft(spectra);

% Perform a very basic phase correction, if specified.
if basicPhaseCorrection
    for voxel = 1:nVoxels
        phase = unwrap(angle(spectra(1,:,voxel)));
        linPhase = (phase(nPoints) - phase(1))/(nPoints-1) * (0:nPoints-1);
        spectra(1,:,voxel) = spectra(1,:,voxel) .* exp(-1i*linPhase);
    end
end

spectra = reshape(spectra, [nPoints spatialSizes nRepetitions]);

% Transpose first two spatial dimensions, if necessary.
if any(Reco.RECO_transposition > 0)
    if sum(Reco.RECO_transposition > 0) == length(Reco.RECO_transposition)
        % entries indicate transposition for all slices
        spectra = permute(spectra, [1 3 2 4 5]);
    else
        % some slices need to be transposed, others not. The result of this
        % operation cannot be stored in a common matrix, because matrix 
        % sizes can be different per slice.
        warning('Transposition could not be applied.');
    end
end

% Calculate abscissa values
ppmAtBaseFreq = rawObj.Method.PVM_FrqWorkPpm(1); % workingChemicalShift
f = transpose(-nPoints/2:(nPoints/2-1)); % column vector, like the spectra

if isEPSI
    abscissaValuesHz = rawObj.Method.SpecBand * f / nPoints;
    abscissaValuesPpm = rawObj.Method.SpecBandPpm * f / nPoints + ppmAtBaseFreq;
else
    abscissaValuesHz = rawObj.Method.PVM_SpecSWH(1) * f / nPoints;
    abscissaValuesPpm = rawObj.Method.PVM_SpecSW(1) * f / nPoints + ppmAtBaseFreq;
end