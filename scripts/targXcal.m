% targXcal: xcalibur TMT targeting list script for generating inclusion 
% lists from MaxQuant evidence.txt files. These can be directly loaded into
% xCalibur. Note that this function is designed for TMT data ONLY as it 
% adds TMT reporter masses to the mass when calculating m/z.
%
% example minimal function call:
%       targXcal(data,'filename')
%
% Where:
%       data = the maxquant evidence.txt imported as a structure or table
%       filename = the output .csv file name for the inclusion list
%
% Optional parameters include: PEP (applies PEP filter), Score (applies
% score filter),  Intensity (mimumum intensity filter),Include (includes 
% only matching Accession numbers), Exclude (excludes only matching 
% Accession numbers), and Unique (allows filtering for unique peptides
% (sequence and z, and will select the most intense).
%
% An optional output is the table of filtered hits generated in the process
% of creating the .csv file.
%
% Ed Emmott, Northeastern U. 2019.

function [targets] = targXcal(data,filename,varargin)

% Input parameter validity checks and defaults
p = inputParser;

% input validity checks
validData = @(x) isstruct(x) | istable(x); % data must be a structure or table
validPEP  = @(x) isnumeric(x) & x >= 0 & x <= 1; % PEP must be a number between 0 and 1.

% list defaults
% Optional Positional parameters
dRTwidth    = 3;  % Default: +/- 3 minute retention time window per peptide
dMZdecimals = 5;  % Default: gives ion m/z to 5dp
dRTdecimals = 2;  % Default: gives retention time in minutes to 1dp

% Optional Name-value parameters
dPEP          = 1;  % Default: PEP filter set at 1 i.e. NO filter applied.
dFilterConRev = true; % Default: removes 'Con','Rev'. Only applied to table data
dScore        = 0;  % Default: no filter
dIntensity    = 0;  % Default: no filter
dRTshift      = 0;
dInclude      = {}; % Default: no filter
dExclude      = {}; % Default: no filter
dUnique       = true;   % Default: if multiple entries for the same ion/z state
% exist, keep the one with the highest intensity.

% list input arguments
addRequired (p,'data'        ,validData);
addRequired (p,'filename'    ,@ischar);
addOptional (p,'rtWidth'     ,dRTwidth,@isnumeric);
addOptional (p,'mzDecimals'  ,dMZdecimals,@isnumeric);
addOptional (p,'rtDecimals'  ,dRTdecimals,@isnumeric);
addParameter(p,'FilterConRev',dFilterConRev,@islogical);
addParameter(p,'PEP'         ,dPEP,validPEP);
addParameter(p,'Score'       ,dScore,@isnumeric); % applies a minimum score threshold
addParameter(p,'Intensity'   ,dIntensity,@isnumeric); % applies a minimum Intensity threshold
addParameter(p,'RTshift'     ,dRTshift,@isnumeric); % subtracts this value from all retention times
addParameter(p,'Unique'      ,dUnique,@islogical); % Want to add several options: first, highest intensity)
addParameter(p,'Include'     ,dInclude,@iscell); % Include only peptides matching these accessions
addParameter(p,'Exclude'     ,dExclude,@iscell); % Exclude peptides matching these accessions. If used together
% with 'Include', 'Include' is applied first, then 'Exclude'.

% Parse input arguments
parse(p,data,filename,varargin{:})

%% Additional check section
% filename .csv check
if strcmp(filename(end-3:end) , '.csv') == 0
    filename = [filename,'.csv'];
    fprintf('The filename given was not a .csv file.\n');
    fprintf(['Your filename has been modified from ',filename(1:end-4),' to ',filename]);
    fprintf('To avoid this in future please use a filename ending in .csv');
end

%% Handling different input data formats (table, structure, and columns); 
% Note: if the spelling or capitalisation of the imported columns changes,
% these will require modification in this section.
    % Col Checker function
    isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);

    %Structure (perl-based) format data
    if isstruct(data) == 1
        data2 = table();
        data2.Accession = data.textdata(:,1);
        data2.Sequences = data.textdata(:,2);%
        data2.PEP = data.PEP; %
        data2.z = data.z; %
        data2.Mass = data.Mass; %
        data2.RT = data.RT; %
        %data2.Score = data.Score; % score column not typically imported by
        % perl script
        data2.Intensity = data.Intensity; %
               
    end

    %Table (readtable-based) format data (default for most users)
    if istable(data) == 1
        data2 = table();
        data2.Accession = data.Proteins;
        data2.Sequences = data.ModifiedSequence;%
        data2.PEP = data.PEP; %
        data2.z = data.Charge; %
        data2.Mass = data.Mass; %
        data2.RT = data.RetentionTime; %
        data2.Score = data.Score; %
        data2.Intensity = data.Intensity; %
               
        % If filtering out 'Con' and 'Rev' hits.
        if p.Results.FilterConRev == true
           revHit = ismember(data.Reverse,'+');
           conHit = ismember(data.PotentialContaminant,'+');
           revCon = sum([revHit,conHit],2) >= 1;
           data2 = data2(~revCon,:); % Removes Rev and Con
        end
             
    end
    
    data = data2;

%% Data filtering (as requested or specified by defaults)

    % PEP filter
    data = data(data.PEP <= p.Results.PEP,:);
    
    % Score Filter (Checks first if filter applied, then if filter exists)
    if p.Results.Score > 0
        if isTableCol(data,'Score') == 0
            fprintf(2,'No score filter could be applied as no score column was identified in the input data');
        else
            data = data(data.Score >= p.Results.Score,:);
        end
    end
    
    % Intensity Filter (Checks first if Intensity column exists)
    if isTableCol(data,'Intensity') == 0
        fprintf(2,'No intensity filter could be applied as no Intensity column was identified in the input data');
    else
        data = data(data.Intensity >= p.Results.Intensity,:);
    end
    
    % Include Filter
    if isempty(p.Results.Include) == 0
        data = data(ismember(data.Accession , p.Results.Include),:);
    end
    
    % Exclude Filter
    if isempty(p.Results.Exclude) == 0
        data = data(~ismember(data.Accession , p.Results.Exclude),:);
    end
    

    % Filter on unique
    if p.Results.Unique == true
    % First generate sequence/charge fusion
        for ii = 1:size(data.Sequences , 1)
            data.seqZ{ii,1} = [data.Sequences{ii,1} , num2str(data.z(ii))];
        end

        % filter for unique ions
        [c , ia , ~] = unique(data.seqZ);
        data2 = data(ia,:); % creates a table with the first entry for each unique seqZ combination.
        for ii = 1:numel(c)
           if sum(ismember(c(ii) , data.seqZ)) > 1
              % copy out the matching bit of table
              [~ , idx ] = ismember(c , data.seqZ);
              data3 = data(idx,:); % make table with only the non-unique entries
              [~ , idx2] = max(data3.Intensity); % identify the entry with highest intensity (ignores nan by default)
              data2(ii,:) = data3(idx2,:); % replaces the entry in data2, with the highest intensity entry              
              clear data3 idx idx2
           end
        end
        data = data2;
        clear c ia data2  
    end

    % Create Targeting table
    targets = table();
    targets.sequences = data.Sequences;
    targets.mass      = data.Mass         ;
    targets.z         = data.z            ;
    targets.RT        = data.RT - p.Results.RTshift; % Applies RTshift if used
    targets.int       = data.Intensity  ;

    % identify number of lysine residues for TMT labelling
    for ii = 1:numel(targets.sequences)
        targets.numK(ii)  = numel(strfind(targets.sequences{ii} , 'K'));
    end

    % number of TMT labels per peptide is numK + 1 (+1 is for the N-terminus);
    targets.numTMT = targets.numK + 1;

    % Calculate theoretical m/z from mass, z and numTMT.
    % adds the mass, number of protons corresponding to z, number of TMT tags
    % for the peptide, and then divides everything by the charge.
    % NOTE: if wishing to use this function for non-TMT data: this is this
    % section that will require modification.
    targets.theoMz   = (targets.mass + (targets.z*1.007276) +...
        (targets.numTMT.*229.162932)) ./ targets.z;
    targets.theoMz = round(targets.theoMz , p.Results.mzDecimals);

    
    % Start end end time
    targets.startRT = round(targets.RT - p.Results.rtWidth , p.Results.rtDecimals);
    targets.endRT   = round(targets.RT + p.Results.rtWidth , p.Results.rtDecimals);

    %% xCalibur-formating
    % Make matrix containing the relevent numeric values
    xcal = [targets.theoMz , targets.z , targets.startRT , targets.endRT ];
    
    % Generate cell array containing sequence and intensity to include in
    % the comments box for DO-MS plugin
    for ii = 1:numel(targets.int)
        comment{ii,1} = [targets.sequences{ii},';',num2str(targets.int(ii))];
    end

    % Output format
    formatString = ['%.',num2str(p.Results.mzDecimals),'f,,,,%g,Positive,%.',...
        num2str(p.Results.rtDecimals),'f,%.',num2str(p.Results.rtDecimals),...
        'f,,,,%s,\n'];

    fileID = fopen(filename , 'w');
    %write headers and format
    fprintf(fileID,'%s%s%s%s%s%s%s%s%s%s%s%s\n',...
    'Mass [m/z],','Formula [M],','Formula type,','Species,','CS [z],','Polarity,',...
    'Start [min],','End [min],','(N)CE,','(N)CE type,','MSX ID,','Comment,');

    for ii = 1:size(xcal,1)
        fprintf(fileID, formatString ,xcal(ii,:),comment{ii,1});
    end

    fclose(fileID);
    
    % Print finished
    fprintf(['\nTargeting list successfully written to file: ',filename,'\n']);
    
    % Print warning warning if inclusion list is oversize for xCalibur.
    if ii >= 5001
        fprintf(2,['WARNING: Maximum inclusion list size exceeded. Inclusion list has ',...
            num2str(ii - 5000),' entries over the 5000 entry maximum for xCalibur.\n']);
        fprintf(2,['While the targeting file has been created, it will not be usable ',...
            'in xCalibur until the number of entries is reduced.\n']);
        fprintf(2,['You should either take steps to reduce targeting list size or consider '...
            'use of alternate targeting software e.g. MaxQuant.live\n']);
        fprintf(2,'which permits the use of longer inclusion lists\n');
    end

end