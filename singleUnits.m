function resp_cids = singleUnits()%% Save units for analysis

% to do: save with *_UnitResponses.mat data
% selected responsive units
vibrotac_nosound = [10121 10334 10400 10441 10457 11212 11247 11259 19153 19265 19287 20277 20290 20303 20306]'; % vibrotac responsive, during plug
pressure_nosound = [10400 10441 10457 11212 11257 19265 19287 19296 20277 20303 20306]'; % pressure responsive, during plug
new_cids = [31321 31322 33110 33117 33205 33238 33245 33272]'; %29198
resp_cids = unique([vibrotac_nosound; pressure_nosound; new_cids]); % vibrotac and/or pressure responsive, exclusively with earplug

%% selection all vibrotactile responsive units: vibrotac and/or pressure responsive, not exclusively with earplug
%selected_cids = [10121 10328 10330 10334 10382 10386 10400 10421 10423 10441 10457 11210 11212 11217 11247 11257 11259 11263 1939 1990 19153 19265 19287 19296 2046 20225 20276 20277 20287 20290 20303 20306]'; 

% selection sound vibrotac responsive units
%nonresp_cids = selected_cids(~ismember(selected_cids, resp_cids));
end