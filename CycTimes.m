function [CycT] = CycTimes(SpkT,StimDur, Mf,SkipVal, SkipMethod)
%CYCTIMES Converts stimulus timing to cycle times for AM stimulus
%   INPUTS:
%   SpkT    = (cell?) array of spike times relative to sound onset
%   StimDur = scalar; duration of stimulus in seconds;
%   Mf      = scalar; Modulation frequency in Hz;
%   SkipVal = scalar; duration in seconds or number of cycles to skip
%   SkipMethod = string; specify how to treat onset response. either 'time'
%   or 'cycles' (default)
% 
%   OUTPUTS:
%   CycT    = cell array of spike times relative to cycles in cycles
%   [0,1)
%

    % Checking values of input
    if ~iscell(SpkT); SpkT = {SpkT};end
    NStim = length(SpkT);
    if length(StimDur) == 1; StimDur = StimDur.*ones(size(SpkT));end
    if length(StimDur) ~= NStim; error('Size of StimDur does not match SpkT.'); end
    
    % default values
    if nargin < 5; SkipMethod = 'cycles';end
    if nargin < 4; SkipVal = 0; end
    
    CycT = cell(size(SpkT));
    if Mf == 0; return; end % return empty cells if modulation frequency is 0
    
    T = 1/Mf; % period in seconds;
    
    switch SkipMethod
        case 'cycles'
            SkipDur = SkipVal * T;
        case 'time'
            SkipDur = SkipVal;
        otherwise
            SkipMethod = 'cycles';
            SkipDur =  SkipVal * T;
    end
    
    for i = 1: NStim
        tSpk = SpkT{i};
        tSpk = tSpk(tSpk > SkipDur & tSpk < StimDur(i)); % select spks during stimulus, skipping onset
        if isempty(tSpk)
            CycT{i} = [];
        else
            CycT{i} = mod(tSpk,T)./T;
        end
        
    end
end

