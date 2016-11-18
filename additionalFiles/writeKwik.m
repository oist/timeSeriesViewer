function writeKwik( filename, data, triggers, triggersChNum )
%WRITEKWIK Writes two kwik files containing data and triggers
%   This function creates two files following the KWIK specifications
%   (here: https://github.com/klusta-team/kwiklib/wiki/Kwik-format#kwik)
%   
%   The file filename.kwd contains the raw data, while the file
%   filename.kwe is the ne containing the trigger information. Two files
%   are used following the variations of the original format specified
%   here: https://open-ephys.atlassian.net/wiki/display/OEW/KWIK+format
%
%   INPUT: 
%     filename the desired filename 
%     data the raw data, passed as a 3D matrix of the format [nChannels x nTrials x nSamples]
%     triggers the triggers information, passed as cell array
%     triggersChNum channel number information for the triggers
%   OUTPUT: none. 2 files will be written on disk
%
%   Author: Stefano.Masneri@brain.mpg.de
%   Date: 18.11.2016

dataFile = [filename '.kwd'];
trigFile = [filename '.kwe'];

% first write the .kwd file
numCh = size(data, 1);
for k = 1:numCh
  h5create(dataFile, ['/recordings/' num2str(k) '/data'], size(data(k, :, :)), ...
    'Datatype', class(data));
  h5write( dataFile, ['/recordings/' num2str(k) '/data'], data(k, :, :));
end

h5writeatt(dataFile,'/','kwik_version',2);

% now for the events
lengthsTriggers = cellfun('length', triggers); % get the length of each element of the cell
trigMat = cell2mat(triggers);                  % concatenate all the cell elements
h5create(trigFile, '/event_types/TTL/events/timesamples', size(trigMat), ...
  'Datatype', class(trigMat));
h5write(trigFile, '/event_types/TTL/events/timesamples', trigMat);

% set the on off value for the triggers
triggersOnOff = zeros(1, length(trigMat));     
triggersOnOff(1:lengthsTriggers(1:2:end)) = 1;
triggersOnOff = logical(triggersOnOff);
h5create(trigFile, '/event_types/TTL/events/userdata/eventID', size(triggersOnOff), ...
  'Datatype', class(triggersOnOff));
h5write(trigFile, '/event_types/TTL/events/userdata/eventID', triggersOnOff);

% set the channel value for the triggers
channelInfo = [];
for k = 1:length(triggersChNum)
  channelInfo(end+1 : end+1+lengthsTriggers(k)) = triggersChNum(k);
end
h5create(trigFile, '/event_types/TTL/events/userdata/event_channels', size(channelInfo), ...
  'Datatype', class(channelInfo));
h5write(trigFile, '/event_types/TTL/events/userdata/eventID', channelInfo);

end