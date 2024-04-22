function out = tetraConverter(input, opts)
% TETRACONVERTER Convert between MNI coordinates on the scalp/brain and
% Tetra Code. Also recognizes and supports conversion of CPC coordinates,
% Beam-style X/Y percentages, and 10-10 EEG electrode locations.
%   out = TETRACONVERTER(input) tries to automatically detect the input
%   type and returns the corresponding scalp MNI coordinates (if input is
%   Tetra Code) or Tetra Code (if input is anything else). Each row in the
%   input would have a corresponding output
%
%   out = TETRACONVERTER(input, 'intype', inputType) allows additional 
%   specification of the input as either 'scalp', 'brain', 'tetra', 'CPC',
%   'XY', or 'EEG'

%   out = TETRACONVERTER(input, 'outtype', outputType) specifies calls for
%   the function to output some combination of 'scalp', 'brain', 'tetra', 
%   'CPC', XY', 'EEG', or 'all' of them. If more than one outputType is 
%   called, the function outputs them as a table
arguments
    input
    opts.intype string = 'auto'
    opts.outtype string = 'auto'
end
assert(any(strcmp(opts.intype, {'auto', 'scalp', 'brain', 'tetra', 'CPC', 'XY', 'EEG'})), ...
    'intype must be ''auto'' (default), ''scalp'', ''brain'', ''tetra'', ''CPC'', ''XY'', or ''EEG''')
assert(all(contains(opts.outtype, {'auto', 'scalp', 'brain', 'tetra', 'CPC', 'XY', 'EEG', 'all'})), ...
    'outtype must be ''auto'' (default), ''scalp'', ''brain'', ''tetra'', ''CPC'', ''XY'', ''EEG'', or ''all''')
if ~isnumeric(input) && ~(ischar(input) || isstring(input))
    error('Input must be a numeric, string, or character array')
end
if strcmp(opts.intype, 'auto')
    if isnumeric(input) % numeric input has to be CPC or MNI
        if width(input)==2 % assume CPC input
            if any(input(:, 1)>=1.1084|input(:, 2)>1)
                opts.intype = 'XY';
                disp('Assuming input is Beam-style X/Y')
            else
                opts.intype = 'CPC';
                disp('Assuming input is CPC')
            end
        elseif width(input)~=3 % otherwise assume MNI input; if not:
            error('Numeric input must be in a Nx3 array (MNI) or Nx2 array (CPC)')
        end
    else
        disp('Assuming input is Tetra Code or EEG')
        opts.intype = 'tetra'; % EEG and tetra are parsed the same way so it doesn't matter here
    end
end

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ref = readtable('tet_code2MNI_lookup_extended.xlsx', 'Sheet','Reference');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

istetra = zeros(1, height(input));
switch opts.intype
    case 'auto'
        disp('Assuming input is MNI scalp or brain')
        assert(width(input)==3, 'MNI input must be in a Nx3 array')
        matchID = findMinDisID([ref{:, {'ScalpX', 'ScalpY', 'ScalpZ'}}; ref{:, {'BrainX', 'BrainY', 'BrainZ'}}], input);
        matchID = mod(matchID, height(ref));
        matchID(~matchID) = height(ref);
    case 'brain'
        assert(width(input)==3, 'MNI input must be in a Nx3 array')
        matchID = findMinDisID(ref{:, {'BrainX', 'BrainY', 'BrainZ'}}, input);
    case 'scalp'
        assert(width(input)==3, 'MNI input must be in a Nx3 array')
        matchID = findMinDisID(ref{:, {'ScalpX', 'ScalpY', 'ScalpZ'}}, input);
    case 'CPC'
        assert(width(input)==2, 'CPC input must be in a Nx2 array')
        assert(all(input>=0, 'all') && all(input<1.1084, 'all'), 'CPC input out of range [0, 1.1084)')
        matchID = findMinDisID(ref{:, {'Pnz', 'Pal'}}, input);
    case 'XY'
        assert(width(input)==2, 'Beam-style X/Y input must be in a Nx2 array')
        assert(all(input>=0, 'all') && all(input<152.0659, 'all'), 'Beam-style X/Y input out of range [0, 152.0659)')
        matchID = findMinDisID(ref{:, {'X_', 'Y_'}}, input);
    case {'EEG', 'tetra'}
        input = string(input);
        assert(width(input)==1, 'EEG or Tetra Code input must be a single column')
        matchID = zeros(1, height(input));
        for row=1:height(input)
            try
                [matchID(row), istetra(row)] = find(strcmpi(input(row), ref{:, {'TetraCode', 'EEG'}}), 1);
            catch
                error(['Match not found for "' char(input(row)) '", check spelling of input'])
            end
        end
end

if strcmp(opts.outtype, 'auto')
    if any(istetra==1)
        opts.outtype = 'scalp';
    else
        opts.outtype = 'tetra';
    end
end

out.scalp = ref{matchID,  {'ScalpX', 'ScalpY', 'ScalpZ'}};
out.brain = ref{matchID,  {'BrainX', 'BrainY', 'BrainZ'}};
out.tetra = string(ref.TetraCode(matchID));
out.CPC = ref{matchID,  {'Pnz', 'Pal'}};
out.XY = ref{matchID,  {'X_', 'Y_'}};
out.EEG = string(ref.EEG(matchID));
out = struct2table(out);
if ~contains(opts.outtype, 'all')
    out = out(:, opts.outtype);
    if width(out)==1
        out = out{:, :};
    end
end
end

function matchID = findMinDisID(ref, targets)
    assert(size(ref, 2) == size(targets, 2), 'Input arrays must have the same width');
    matchID = zeros(1, height(targets));
    for row = 1:size(targets, 1)
        [~, matchID(row)] = min(sum((ref - targets(row, :)).^2, 2));
    end
end