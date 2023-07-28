function [ structOut ] = bruker_combineStructs( structA, structB )
%COMBINESTRUCTS Summary of this function goes here
%   combines structA and structB to a struct: structOut

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

names=fieldnames(structB);
for i=1:length(names)
    if ~isfield(structA, names{i})
        structA.(names{i})=structB.(names{i});
    else
        disp(['Warning: same field ', names{i} ' in both structs. No overwrite !']);
    end

end
structOut=structA;
