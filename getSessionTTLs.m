function [sessions_TTLs, sessions_TTLs_details] = getSessionTTLs(messagesPath, rec_samples, Fs)
% sessions_TTLs_details contains cell array: text from message center,
% start(1)/end(0) code, corresponding sample nr, time since start recording
% sessions_TTLs contains array with: event nr, start(1)/end(0) code, corresponding sample nr

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

% add sample nr of corresponding message
ttl_samples = [num2cell(message_samples(1:size(msgs, 1)))]; % sample nr
sessions_TTLs_details = [msgs, ttl_samples];

ttls_norm = num2cell((message_samples - rec_samples(1))/Fs);
sessions_TTLs_details = [sessions_TTLs_details, ttls_norm]; % time since recording onset

%filename = sprintf('M%.2i_S%02d-%02d_OE_TTLs', str2double(stimuli_parameters.Par.MouseNum), relevant_sessions(1), relevant_sessions(2));
%save(fullfile(OutPath, filename), "sessions_TTLs")


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