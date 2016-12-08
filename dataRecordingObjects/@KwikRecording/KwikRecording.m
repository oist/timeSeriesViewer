classdef KwikRecording < dataRecording
  %KWIKRECORDING Reads a recording in the Kwik format
  % This class is used to read kwd and kwe files containing raw data and
  % trigger information from a kwik recording. The class DOES NOT implement
  % the whole kwik data specification but it only focuses on data and
  % triggers at the moment.
  %
  % AUTHOR: stefano.masneri@brain.mpg.de
  % DATE: 21.11.2016
  % 
  % TODO:
  %   - check how to deal with channel info for triggers
  %   - add functionalities
  
  properties
    numRecordings;   % number of recordings in a single .kwd file
    timestamps;      % timestamps information for each channel
    triggerFilename; % name of the *.kwe file containing trigger info
    info;            % additional information provided in the file
    bitDepth;        % number of bits used to store data
    datatype;        % class of data in the recording
    sample_ms;
    fullFilename;    % path + name
    recordingNames;  % names of all recordings
  end
  
  properties (Constant = true)
    pathToData = '/recordings/';
    
    pathToTriggerData = '/event_types/TTL/events/time_samples'; % where the triggers are stored in the .kwe file
    pathToTriggerOnOff = '/event_types/TTL/events/user_data/eventID';
    pathToTriggerChannel = '/event_types/TTL/events/user_data/event_channels';
        
    %must define them, bc abstract in base class
    defaultLocalDir='C:\Users\Tulip\Documents\Academic\Post-Doc\Experiments'; %Default directory from which search starts
    signalBits = 16; %the quantization of the sampling card
    numberOfCharFromEndToBaseName=7;
  end
  
  methods
    
    function [V_uV ,t_ms]=getData(obj, channels, startTime_ms, window_ms)
      %GETDATA Extract Kwik recording data from file to memory
      %Usage: [V_uV,t_ms] = obj.getData(channels,startTime_ms,window_ms);
      %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
      %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
      %        window_ms - a scalar [1x1] with the window duration [ms].
      %Output: V_uv - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
      %        t_ms - A time vector relative to recording start (t=0 at start)
      
      windowSamples = double(round(double(window_ms) / obj.sample_ms(1)));
      nWindows = numel(startTime_ms);
      startTime_ms = round(startTime_ms/obj.sample_ms(1))*obj.sample_ms(1);
      startElement = double(round(startTime_ms/obj.sample_ms(1)));
      if 0 == startElement || Inf == startElement
        startElement = 1;
      end
      %window_ms = windowSamples * obj.sample_ms(1);
      
      if isempty(channels) %if no channels are entered, get all channels
        channels=obj.channelNumbers;
      end
      nCh = numel(channels);
      
      V_uV = zeros(nCh, nWindows, windowSamples, obj.datatype); %initialize waveform matrix
      
      for k = 1:length(channels)
        for m = 1:numel(startElement)
          V_uV(channels(k), :, :) = h5read(obj.fullFilename, [obj.recordingNames{1} '/data'], ...
            [k startElement(m)], [1 windowSamples]);
%         try
%           V_uV(channels(k), :, :) = h5read(obj.fullFilename, [obj.channelNames{k} '/data'], ...
%             [1 1 startElement], [1 nWindows windowSamples]);
%         catch
%           V_uV(channels(k), 1, :) = h5read(obj.fullFilename, [obj.channelNames{k} '/data'], ...
%             [1 startElement], [1 windowSamples]);
        end
      end
      
      if obj.convertData2Double
          V_uV=double(V_uV);
          for k = 1:size(V_uV, 1)
              V_uV(k, :, :) = V_uV(k, :, :) * obj.MicrovoltsPerAD(k);
          end
      end
      
      if nargout==2
        t_ms=(1:windowSamples)*(1e3/obj.samplingFrequency(1));
      end

    end
    
    function [T_ms, chNumber] = getTrigger(obj, ~, ~, ~)
      %GETTRIGGER Extract triggers from file Kwik recording
      % Read Data from *.kwe files (NeuroLynx) - files containing the data
      %Usage : [T_ms] = obj.getTrigger(name,startTime_ms,endTime_ms)
      %Input : name - which bit to extract for time stamps (out of 8,default = first bit, 1)
      %        startTime_ms - start time [ms].
      %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
      %
      %Output: T_ms - trigger times [ms] - different triggers are arranged in a cell array
      
      T_ms = cell(1, 2);
      allTriggers = h5read(obj.triggerFilename, obj.pathToTriggerData);
      onOff =  h5read(obj.triggerFilename, obj.pathToTriggerOnOff);
      T_ms{1} = allTriggers(onOff == 1);
      T_ms{2} = allTriggers(onOff == 0);
      chNumber = h5read(obj.triggerFilename, obj.pathToTriggerChannel);
    end
  end
  
  methods (Hidden = true)
    
    %class constructor
    function obj = KwikRecording(recordingFile)
      %get data files
      if nargin == 0
        recordingFile=[];
      elseif nargin>1
        disp('KwikRecording: Object was not constructed since too many parameters were given at construction');
        return;
      end
      
      obj = obj.getRecordingFiles(recordingFile, 'kwd');
      
      % Find the .kwe file
      %[~, name, ~] = fileparts(obj.recordingName);
      obj.fullFilename = fullfile(obj.recordingDir, obj.recordingName);
      triggerFile = dir([obj.recordingDir filesep '*.kwe']);
      if isempty(triggerFile)
        error('KwikRecording: Cannot file .kwe file')
      elseif length(triggerFile) > 1
        warning('KwikRecording: Multiple .kwe file found! using the first one')
      end
      obj.triggerFilename = fullfile(obj.recordingDir, triggerFile.name);
      
      fileInfo = h5info(obj.fullFilename, obj.pathToData);
      
      obj.numRecordings = length(fileInfo.Groups);
      if obj.numRecordings > 1
        warning('KwikRecording: file contains multiple recordings.')
      end
      obj.recordingNames = cell(1, obj.numRecordings);
      for k = 1:obj.numRecordings
        obj.recordingNames{k} = fileInfo.Groups(k).Name;
      end
      
      % Now extract metadata
      try
        obj.MicrovoltsPerAD = double(h5readatt(obj.fullFilename, [obj.recordingNames{1} '/application_data'], 'channel_bit_volts'));
      catch
        obj.MicrovoltsPerAD = double(h5read(obj.fullFilename, [obj.recordingNames{1} '/application_data/channel_bit_volts']));
      end
      
      obj.channelNumbers = 1:length(obj.MicrovoltsPerAD);
      obj.channelNames = cellfun(@(x) num2str(x),mat2cell(obj.channelNumbers,1,ones(1,numel(obj.channelNumbers))),'UniformOutput',0);
      
      try
        obj.samplingFrequency = double(h5readatt(obj.fullFilename, [obj.recordingNames{1} '/application_data'], 'channel_sample_rates'));
      catch
        obj.samplingFrequency = double(h5read(obj.fullFilename, [obj.recordingNames{1} '/application_data/channel_sample_rates']));
      end
      
      try
        obj.timestamps = double(h5read(obj.fullFilename, [obj.recordingNames{1} '/application_data/timestamps']));
      catch
        disp('KwikRecording: timestamps information not available')
      end

      obj.bitDepth = double(h5readatt(obj.fullFilename, obj.recordingNames{1}, 'bit_depth'));

      try
        obj.recordingDuration_ms = double(h5readatt(obj.fullFilename,[obj.recordingNames{1} '/data'], 'valid_samples'));
        obj.recordingDuration_ms = 1000 * obj.recordingDuration_ms ./ obj.samplingFrequency;
        obj.recordingDuration_ms = max(obj.recordingDuration_ms);
      catch
        try
          obj.recordingDuration_ms = double(h5readatt(obj.fullFilename, '/', 'recordingDuration'));
        catch
          obj.recordingDuration_ms = 0;
        end
      end

      obj.sample_ms = 1e3 / obj.samplingFrequency(1);
      
      try
        obj.datatype = h5readatt(obj.fullFilename, '/', 'datatype');
      catch
        obj.datatype = ['int' num2str(obj.bitDepth)];
      end
      
    end
    
    function delete(obj) %do nothing
    end
    
  end
  
end

