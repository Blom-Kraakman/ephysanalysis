function [sessions_TTLs, sessions_TTLs_details] = getSessionTTLs(messagesPath, rec_samples, Fs)
% sessions_TTLs_details contains cell array: text from message center, 
% start(1)/end(0) code, corresponding sample nr, time since start recording
% sessions_TTLs contains array with: event nr, start(1)/end(0) code, corresponding sample nr

% messagesPath = '\\store\department\neuw\shared\Aaron Wong\Data\EphysRecordings\M7\M07_2024-02-01_14-09-33\Record Node 103\experiment1\recording1\events\MessageCenter\';
% recPath = '\\store\department\neuw\shared\Aaron Wong\Data\EphysRecordings\M7\M07_2024-02-01_14-09-33\Record Node 103\experiment1\recording1\continuous\Intan-100.Rhythm Data-A\';

% open files
message_text = [messagesPath 'text.npy']; % session TTLs
message_samples = readNPY([messagesPath 'sample_numbers.npy']); % session TTLs

% open text messages
msgs = getMessageText(message_text);

% make (cell)array to contain session TTL variables
sessions_TTLs = (1:size(msgs, 1))';
for i = 1:size(msgs, 1)
    if contains(msgs{i,1}, 'start')
        msgs{i,2} = 1;
        sessions_TTLs(i,2) = 1;
    elseif contains(msgs{i,1}, 'end')
        msgs{i,2} = 0;
        sessions_TTLs(i,2) = 0;
    else
        error('Error in getMessageText, check file')
    end
end
sessions_TTLs(:,3) = message_samples(1:length(sessions_TTLs)); % sample nr

% add sample nr of corresponding message
ttl_samples = [num2cell(message_samples(1:size(msgs, 1)))]; % sample nr
sessions_TTLs_details = [msgs, ttl_samples];

ttls_norm = num2cell((message_samples - rec_samples(1))/Fs);
sessions_TTLs_details = [sessions_TTLs_details, ttls_norm]; % time since recording onset

function msgs = getMessageText(message_text)
% open messages file
fid = fopen(message_text, 'r');
msg = textscan(fid, '%s', 'delimiter', {'\n'});
fclose(fid);

% split messages by '0' chars (these are null)
msgs = strsplit(msg{1}{2}, '\0')';

% remove empties
msgs = msgs(~cellfun('isempty', msgs));
end
end