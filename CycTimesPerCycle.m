function [CycT,edges] = CycTimesPerCycle(SpkT,StimDur, Freq, Delay)
%CYCTIMES Converts stimulus timing to cycle times for cyclical stimulus
%   INPUTS:
%   SpkT    = 1-D (cell) array of length NStim of spike times relative to stim onset
%   StimDur = scalar; duration of stimulus in seconds;
%   Freq    = scalar; Modulation frequency in Hz;
%   Delay   = scalar; duration in seconds to shift the onset of the first
%   cycle
% 
%   OUTPUTS:
%   CycT    = (NStim x NCyc) cell array of spike times 
%

    % Checking values of input
    if ~iscell(SpkT); SpkT = {SpkT};end
    NStim = length(SpkT);
    if length(StimDur) == 1; StimDur = StimDur.*ones(size(SpkT));end
    if length(StimDur) ~= NStim; error('Size of StimDur does not match SpkT.'); end
    
    % default values
    if nargin < 5; SkipMethod = 'cycles';end
    if nargin < 4; SkipVal = 0; end
    
    if Freq == 0 
        T = StimDur; 
    else
        T = 1/Freq; % period in seconds;
    end % return empty cells if modulation frequency is 0
    NCyc = StimDur./ T;
    CycT = cell(NStim,max(NCyc));
    
    edges = Delay:T:(Delay+StimDur);

    for stim_idx = 1: NStim
        tSpk = SpkT{stim_idx};
        Y = discretize(tSpk,edges);
        if any(~isnan(Y))
            for cyc_idx = 1:max(Y)
                CycT{stim_idx,cyc_idx} = tSpk(Y==cyc_idx);
            end
        end
    end
end

