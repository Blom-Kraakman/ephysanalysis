function [sessions_TTLs, sessions_TTLs_variables] = getSessionTTLs(messagesPath, rec_samples, Fs, skip_sessions)
% sessions_TTLs contains array with: event nr, start(1)/end(0) code,
% corresponding sample nr, time since start recording
% sessions_TTLs_variables contains cell array of labels

% open files
message_text = [messagesPath 'text.npy']; % session TTLs text
message_samples = readNPY([messagesPath 'sample_numbers.npy']); % session TTLs sample

% open text messages
msgs = getMessageText(message_text);

% make (cell)array to contain session TTL variables
sessions_TTLs = (1:size(msgs, 1))';
for i = 1:size(msgs, 1)
    if contains(msgs{i,1}, 'start')
        msgs{i,2} = 1;
        sessions_TTLs(i,2) = 1;
        sessions_TTLs(i,1) = str2double(extractAfter(msgs(i), 'start s'));
    elseif contains(msgs{i,1}, 'end')
        msgs{i,2} = 0;
        sessions_TTLs(i,2) = 0;
        sessions_TTLs(i,1) = str2double(extractAfter(msgs(i), 'end s'));
    else
        error('Error in getMessageText, check file')
    end

end
sessions_TTLs(:,3) = message_samples(1:length(sessions_TTLs)); % sample nr
ttls_norm = double(message_samples - rec_samples(1))/Fs;
sessions_TTLs(:,4) = ttls_norm;

% keep only sessions in recording stretch
if ~isempty(skip_sessions)
    idx_remove = find(sum(sessions_TTLs(:,1) == skip_sessions, 2));
    sessions_TTLs(idx_remove, :) = [];
end

sessions_TTLs_variables = {'session nr', 'TTL on/off', 'TTL (sample nr)', 'TTL (sec since start)'};

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