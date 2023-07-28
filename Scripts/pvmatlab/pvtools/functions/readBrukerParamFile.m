function [paramStruct,headers]=readBrukerParamFile(filename,updateWithOutputFiles)

% Reads Bruker JCAMP parameter files.
% In particular, this function is used to read acqp, method, reco and 
% visu_pars files.
%
% Usage: [paramStruct,headers]=readBrukerParamFile(filename,updateWithOutputFiles)
%
% paramStruct           : Structure containing parameter values.
%                         The parameter names are derived from the JCAMP 
%                         tags.
% headers               : Cell array of strings containing the parsed lines 
%                         from the file header.
% filename              : Name of the parameter file to be read (full path).
% updateWithOutputFiles : [optional] If false, only read the original 
%                         parameter file as specified, do not try to update 
%                         parameters by interpreting any accompanying .out 
%                         result file. The default is true.
%
% Please note: Starting from ParaVision 360 V2.1, all parameter files may 
% have an accompanying result / output parameter file. For example, next to
% the acqp file, a file called acqp.out may be present. Parameters from the
% result / output parameter files always take precedence over the values 
% stored in the original parameter file. The readBrukerParamFile function 
% implicitly reads the result / output file (if present) and updates the 
% parameters accordingly. Set the updateWithOutputFiles argument to false 
% to skip reading the output file. To switch off output files globally, set 
% updateWithOutputFiles=false in the base workspace.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2013-2021
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('updateWithOutputFiles', 'var'))
    try
        updateWithOutputFiles = evalin('base','updateWithOutputFiles');
    catch
        updateWithOutputFiles = true;
    end
end

% Initialize empty parameter structure.
paramStruct = struct;

% Read the original parameter file, fill in paramStruct and headers.
[paramStruct,headers]=readBrukerParamFileUpdate(filename,paramStruct);

outputFileName = [filename '.out'];

if (updateWithOutputFiles && exist(outputFileName, 'file'))
    % Read the result / output parameter file and update the paramStruct.
    [paramStruct]=readBrukerParamFileUpdate(outputFileName,paramStruct);
end


%% readBrukerParamFileUpdate - read a parameter file and update or extend a parameter struct.
% -----------------------------------------------------------------

function [paramStruct,headers]=readBrukerParamFileUpdate(filename,paramStructIn)

% Reads Bruker JCAMP parameter files, updating an existing parameter struct.
%
% Usage: [paramStruct,headers]=readBrukerParamFileUpdate(filename,paramStructIn)
%
% paramStruct   : Structure containing parameter values from paramStructIn, 
%                 updated from the parameter file or extended by additional
%                 parameters found in the parameter file.
%                 The parameter names are derived from the JCAMP tags.
% headers       : Cell array of strings containing the parsed lines from
%                 the file header.
% filename      : Name of the parameter file to be read (including path).
% paramStructIn : Structure containing parameter values to be updated or
%                 extended.

%% Opening and first tests


% Open parameter file
try
    fid = fopen(filename,'r');
catch
    fid = -1;
end
if fid == -1
    error('Cannot open parameter file. Problem opening file %s.',filename);
end

% Init:
count = 0;
line=fgetl(fid);
headers=cell(15,2);

% generate HeaderInformation:
while ~(strncmp(line, '##$', 3))
    count = count+1;
    
    % Retreive the Labeled Data Record    
    [field,rem]=strtok(line,'=');
    
    if ~(strncmp(line, '##', 2)) % it's a comment
        headers{count,2} = strtrim(strtok(line,'$$'));
        if strncmp(headers{count,2}, filesep, 1)
            headers{count,1}='Path';
        elseif strncmp(headers{count,2}, 'process', 7)
            headers{count,1}='Process';
            headers{count,2} = headers{count,2}(9:end);
        else
            pos=strfind(headers{count,2}(1:10), '-');
            if strncmp(headers{count,2},'Mon',3)||strncmp(headers{count,2},'Tue',3)||strncmp(headers{count,2},'Wed',3)||strncmp(headers{count,2},'Thu',3)|| ...
                strncmp(headers{count,2},'Fri',3)||strncmp(headers{count,2},'Sat',3)|| strncmp(headers{count,2},'Sun',3)|| ...
                ( strncmp(headers{count,2}, '20', 2) && length(pos)==2 )
                    headers{count,1}='Date';
            end

        end
    else % it's a variable with ##    
        % Remove whitespaces and comments from  the value
         value = strtok(rem,'=');
         headers{count,1}=strtrim(strtok(field,'##'));
         headers{count,2} = strtrim(strtok(value,'$')); % save value without $ 
    end
        
    
    line=fgetl(fid);
end
headers=headers(1:count,:);

% Check if using a supported version of JCAMP file format
clear pos;
pos=find(strcmpi(headers(:,1), 'JCAMPDX')==1);
if isempty(pos)
    pos=find(strcmpi(headers(:,1), 'JCAMP-DX')==1);
end
if ~isempty(pos) && length(pos)==1
    version = sscanf(headers{pos,2},'%f');
    if (version ~= 5)&&(version ~= 4.24)
        warning(['JCAMP version %f is not supported. '...
                 'The function may not behave as expected.'],version);
    end 
else
    error('Your fileheader is not correct')
end


%% Reading in parameters

% Initialization of parameter struct

paramStruct=paramStructIn;

% set start bool, because line is already read
first_round=true;

% Keep track of "END" tag to prevent reading incomplete files.
% However, allow reading further than "END", since some customers 
% simply concatenate several parameter files and expect that to work ...
sawEND = 0;

% Loop for reading parameters
while ~feof(fid)
    
    % Reading in line
    if ~first_round
        line = getnext(fid, sawEND);
    end
    first_round=false;

    try
        [cell_field]=textscan(line,'##%s %s','delimiter','=');
        field=cell_field{1};
        value=cell_field{2};
    catch
        continue;
    end
    
    % Checking if field present and removing proprietary tag
    try
        field = field{1};
        if strncmp(field,'$',1)
            field=field(2:end);
        end
    catch
        continue;
    end

    % Checking if value present otherwise value is set to empty string
    try
        value = value{1};
    catch
        value = '';
    end
    
    % Checking for END tag
    if strcmp(field,'END')
        sawEND = 1;
        continue;
    end
    
    % Checking if value is a struct, an array or a single value
    if strncmp(value,'( ',2)
        if(strncmp(value, '( <',3)||strncmp(value,'(<',2))
            %is it an dynamic enum ?
            if(strncmp(value, '( <',3))
                value=getDynEnumValues(fid,[2],value);
            elseif (strncmp(value,'(<',2))
                value=getDynEnumValues(fid,[2],value);
            end
                
        else       
            sizes=textscan(value,'%f','delimiter','(,)');
            sizes=sizes{:};
            sizes=sizes(2:end).';
            value=getArrayValues(fid,sizes,'');
            try
                if ~ischar(value{1}) || length(value)==1% possible datashredding with e.g. {'string1', 'string2'}
                    value=cell2mat(value);
                end
            end
        end
        
    elseif strncmp(value,'(',1)
        value=getArrayValues(fid,1,value);

    else
        testNum = str2double(value);
        if ~isnan(testNum)
            value=testNum;
        else
            testString = parseStringOrStringArray(value);
            if ~isempty(testString)
                value=testString;
            end
        end
    end
    
    % Generating struct member
    paramStruct.(field)=value;
    
end

fclose(fid);




%% getnext - gets the next valid line from the file,
% ignores comments and custom fields.
% -----------------------------------------------------------------
function data = getnext(fid, sawEND)
% data   : line data as a string
% fid    : file identifier
% sawEnd : was the END token seen yet?

% Read line
data = fgetl(fid);

% Throwing away comments and empty lines
while strncmp(data,'$$',2)||isempty(data)
    data = fgetl(fid);
end
% Checking for unexpected end of file
if data<0 
    if (sawEND == 0)
        error('Unexpected end of file: Missing END Statement')
    end
else
    data=commentcheck(data);
end




%% getArrayValues - reads an array of values from the file.
% -----------------------------------------------------------------
function values=getArrayValues(fid,sizes,totalData)
% values    : array values read from file
% fid       : file identifier
% sizes     : expected sizes of the array
% totalData : array data already read in 

% Read until next JCAMP tag, comment or empty line; error if unexpected end
% of file occurs
pos = ftell(fid);
data = fgetl(fid);
inString = 0;
while ~(strncmp(data,'##',2)||strncmp(data,'$$',2)||isempty(data))
    thisLine = replace(data, '\<', '');
    thisLine = replace(thisLine, '\>', '');
    inString = inString + length(strfind(thisLine, '<')) - length(strfind(thisLine, '>'));
    %special-case: \ at end of line -> should be \n
    if(strcmp(data(end),'\'))
        totalData=[totalData data 'n '];
    else
        data=commentcheck(data);
        totalData=[totalData data];
        if inString == 0
            % only insert spaces at line breaks outside strings
            totalData = [totalData ' '];
        end
    end
    pos = ftell(fid);
    data = fgetl(fid);
end
fseek(fid, pos, 'bof');
if data<0
    error('Unexpected end of file: Missing END Statement')
end

% Removing whitespaces at the edge of strings
totalData=strrep(totalData,'< ','<');
totalData=strrep(totalData,' >','>');

% Unpack compressed values. For example, replace @4*(0) with 0 0 0 0
expression = '@(\d+?)\*\((.+?)\)';
replacement = '${repmat([$2 '' ''],1,str2num($1))}';
totalData = regexprep(totalData,expression,replacement);
        
% Checking if array is a string or string array ...
if strncmp(totalData,'<',1)
    values=parseStringOrStringArray(totalData);

% ... or an array of structs ...
elseif strncmp(totalData,'(',1)
    values=parseStructArray(totalData);
    
% ... or a simple array (most frequently numeric)
else
    
    values=textscan(totalData,'%s');
    totalStatus=true;
    for count=1:length(values);
        testNum = str2double(values{count});
        if ~isnan(testNum)
            values{count}=testNum;
        else
            totalStatus = false;
        end
    end
    if totalStatus
        values=cell2mat(values);
    end
    try
        values=values{:};
    end

    %flip sizes, since the 'fastest' dimension in memory is the first on in
    %matlab, but the last one in paravision
    values=reshape(values,[fliplr(sizes) 1]);
    
    %flip dimensions back to keep convention of paravision
    values=permute(values,[ndims(values):-1:1]); 

end

%% getDynEnumValues - reads an array of dynamic enums from the file.
% -----------------------------------------------------------------
function values=getDynEnumValues(fid,sizes,totalData)
% values    : array values read from file
% fid       : file identifier
% sizes     : expected sizes of the array
% totalData : array data already read in 

% Read until next JCAMP tag, comment or empty line; error if unexpected end
% of file occurs
pos = ftell(fid);
data = fgetl(fid);
while ~(strncmp(data,'##',2)||strncmp(data,'$$',2)||isempty(data))
    totalData=[totalData data ' '];
    pos = ftell(fid);
    data = fgetl(fid);
end
fseek(fid, pos, 'bof');
if data<0
    error('Unexpected end of file: Missing END Statement')
end

%string shoud be '( <bla> , <blub> )'
not_empty=true;
count=1;
while (not_empty)
    %Remove '( '
    [trash, right]=strtok(totalData,'<');
    right=right(2:end);%remove <
    [left,right]=strtok(right,'>');
    values{count, 1}=left;
    [trash,right]=strtok(right,'<');
    right=right(2:end);%remove <
    [left,right]=strtok(right,'>');
    values{count, 2}=left;
    [trash,totalData]=strtok(right,'<');
    not_empty=~isempty(totalData);
    clear trash;
    
end

%% commentcheck - filter comments
% -----------------------------------------------------------------
function data=commentcheck(data)
%string contains $$?
pos=strfind(data, '$$');
if ~isempty(pos)
    %string contains also < or >?
    if(    ( ~isempty(strfind(data,'<')) ) || ( ~isempty(strfind(data,'>')) )     )
        pos_strstart=strfind(data, '<');
        pos_strend=strfind(data, '>');
        %if $$ is between < > its ok, but if its not, w3e have to remove
        %the rest of line
       stop_check=false; % set true when comment found
        for i=1:length(pos)
            if(~stop_check)
                comment_ok=false;
                for j=1:min([length(pos_strstart), length(pos_strend)])
                    if (pos(i)>pos_strstart(j) && pos(i)<pos_strend(j))
                        comment_ok=true; % set comment_ok to false, when $$ is between of of the <> pairs
                    end
                end
                if ~comment_ok
                    disp(['"', data(pos(i):end), '" removed as comment']);
                    data=data(1:pos(i)-1);
                end
            end
        end
        
    else %string contains only $$ and no <>
        disp([data(pos(1):end), ' removed as comment']);
        data=data(1:pos(1)-1);
        
        
    end
end


%% parseStringOrStringArray - parse string representing a JCAMP string or string array.
% -----------------------------------------------------------------
% Convert the JCAMP-delimited string(s) to a MATLAB string or split into a
% cell array of MATLAB strings.
function outStrings = parseStringOrStringArray(stringToParse)

stringMasked = maskString(stringToParse);
openAng = strfind(stringMasked,'<');
closeAng = strfind(stringMasked,'>');
nStrings = length(openAng);

if nStrings ~= length(closeAng) || sum(openAng < closeAng) ~= nStrings
    error(['Error reading string or string array. ' ...
        'Angle bracket imbalance in: \n', stringToParse]);
end

switch nStrings
    case 0
        % return an empty cell array
        outStrings = {};
    case 1
        % return a string
        outStrings = stringToParse((openAng(1)+1):(closeAng(1)-1));
        % revert escapes for angle brackets needed to store in JCAMP
        outStrings = strrep(outStrings,'\>','>');
        outStrings = strrep(outStrings,'\<','<');
    otherwise
        % return a non-empty cell array of strings
        outStrings = cell(1, nStrings);
        for i = 1 : length(openAng)
            thisString = stringToParse((openAng(i)+1):(closeAng(i)-1));
            thisString = strrep(thisString,'\>','>');
            thisString = strrep(thisString,'\<','<');
            outStrings{i} = thisString;
        end
end

%% parseStructArray - parse string representing an array of structs.
% -----------------------------------------------------------------
% Split the input string into a cell array of struct elements.
function elements = parseStructArray(stringToParse)

maskedString = maskString(stringToParse);

elements = {};
elementIdx = 1;

% Split the input into top-level struct elements, by finding matching outer
% round brackets
openBr = strfind(maskedString,'(');
closeBr = strfind(maskedString,')');
nextOpenBr = 1;
while nextOpenBr <= length(openBr)
    start = openBr(nextOpenBr);
    matchingCloseBr = nextOpenBr;
    % look for closing bracket index "matchingCloseBr", for which there are 
    % "matchingCloseBr" opening brackets with smaller index
    while matchingCloseBr <= length(closeBr) && sum(openBr < closeBr(matchingCloseBr)) ~= matchingCloseBr
        matchingCloseBr = matchingCloseBr + 1;
    end
    if matchingCloseBr > length(closeBr)
        error(['Error reading struct array. Bracket imbalance in: \n', stringToParse]);
    end
    stop = closeBr(matchingCloseBr);

    elementString = stringToParse((start+1):(stop-1));
    elementMasked = maskedString((start+1):(stop-1));    
    element = convertStructStringToCellArray(elementString, elementMasked);

    elements(:,elementIdx) = element;

    nextOpenBr = matchingCloseBr + 1;    
    elementIdx = elementIdx + 1;
end

%% convertStructStringToCellArray - convert a struct to a cell array.
% ------------------------------------------------------------------
% Convert a string representing a struct to a cell array of struct members
% Note that nested members that are themselves structs are flattened.
function element = convertStructStringToCellArray(structString, structStringMasked)

if isempty(strtok(structString))
    % empty or contains only spaces
    element={};
    return;
end

% flatten struct elements by removing remaining round brackets
brackets = [strfind(structStringMasked, '(') strfind(structStringMasked, ')')];
structString(brackets) = ' ';

% find indices of struct member delimiters in masked string
commas = strfind(structStringMasked, ',');
nMembers = length(commas) + 1;
element=cell(nMembers, 1);

for memberIdx = 1:nMembers
    if memberIdx == 1 
        start = 1;
    else
        start = commas(memberIdx-1) + 1;
    end

    if memberIdx == nMembers
        stop = length(structString);
    else
        stop = commas(memberIdx) - 1;
    end

    memberString = structString(start:stop);
    memberStringMasked = structStringMasked(start:stop);

    % Check for string member in masked string
    if contains(memberStringMasked,'<')
        member = parseStringOrStringArray(memberString);
    else
        [member, converted] = str2num(memberString);
        if ~converted
            % could not be converted to numeric. Return string without
            % whitespace
            member = strip(memberString);
        end
    end

    element{memberIdx} = member;
end

%% maskString - mask the content of a JCAMP string.
% ------------------------------------------------------------------
% Create a masked version of the angle-bracket-delimited JCAMP string, 
% in which the content of string(s) is replaced with XXXX. This is to avoid
% that round brackets, escaped angle brackets and commas in the string 
% content do not interfere with later searches for real delimiters.
function maskedString = maskString(jcampString)

maskedString = jcampString;
maskedString = replace(maskedString,'\>','XX');
maskedString = replace(maskedString,'\<','XX');
openAng = strfind(maskedString,'<');
closeAng = strfind(maskedString,'>');
if length(openAng) ~= length(closeAng) || sum(openAng < closeAng) ~= length(openAng)
    error(['Error masking JCAMP string: imbalance in string delimiters:\n', ...
        jcampString]);
end
for i = 1:length(openAng)
    maskedString((openAng(i)+1):(closeAng(i)-1)) = 'X';
end
