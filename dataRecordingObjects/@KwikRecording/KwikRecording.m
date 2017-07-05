classdef KwikRecording < dataRecording
  %KWIKRECORDING Reads a recording in the Kwik format
  % This class is used to read kwd and kwe files containing raw data and
  % trigger information from a kwik recording. The class DOES NOT implement
  % the whole kwik data specification but it only focuses on data and
  % triggers at the moment.
  %
  % AUTHOR: stefano.masneri@brain.mpg.de & Mark Shein-Idelson
  % DATE: 18.04.2017
  % 
  % TODO:
  %   - check how to deal with channel info for triggers
  %   - add functionalities
  
  properties
    numRecordings;   % number of recordings in a single .kwd file
    timestamps;      % timestamps information for each channel
    triggerFilename; % name of the *.kwe file containing trigger info
    bitDepth;        % number of bits used to store data
    sample_ms;
    fullFilename;    % path + name
    recNameHD5;  % names of all recordings
    dataLength;      % total samples in the data
    info;            % information on recording file
    lengthInfo;      % information on recording file length
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
      if startElement == 0 || Inf == startElement
        startElement = 1;
      end
      %window_ms = windowSamples * obj.sample_ms(1);
      
      if isempty(channels) %if no channels are entered, get all channels
        channels=obj.channelNumbers;
      end
      nCh = numel(channels);
      
      V_uV = zeros(nCh, nWindows, windowSamples, obj.datatype); %initialize waveform matrix
      
      % Speed up if all channels are consecutives
      if all(diff(channels)==1)
        for m = 1:numel(startElement)
          if startElement <= -windowSamples
            %do nothing, return all zeros
          elseif startElement(m) < 1
            V_uV(:, :, -startElement(m)+1 : end) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
              [channels(1) 1], [length(channels) windowSamples + startElement(m)]);
          elseif startElement >= obj.dataLength
            startElement = obj.dataLength - windowSamples;
            V_uV(:, :, :) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
              [channels(1) startElement(m)], [length(channels) windowSamples]);
          elseif startElement + windowSamples > obj.dataLength
            V_uV(:, :, 1:obj.dataLength-startElement) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
              [channels(1) startElement(m)], [length(channels) obj.dataLength - startElement]);
          else
            V_uV(:, :, :) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
              [channels(1) startElement(m)], [length(channels) windowSamples]);
          end
        end
      else
        for k = 1:length(channels)
          for m = 1:numel(startElement)
            
            if startElement <= -windowSamples
              %do nothing, return all zeros
            elseif startElement(m) < 1
              V_uV(k, :, -startElement(m)+1 : end) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
                [channels(k) 1], [1 windowSamples + startElement(m)]);
            elseif startElement >= obj.dataLength
              startElement = obj.dataLength - windowSamples;
              V_uV(k, :, :) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
                [channels(k) startElement(m)], [1 windowSamples]);
            elseif startElement + windowSamples > obj.dataLength
              V_uV(k, :, 1:obj.dataLength-startElement) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
                [channels(k) startElement(m)], [1 obj.dataLength - startElement]);
            else
              V_uV(k, :, :) = h5read(obj.fullFilename, [obj.recNameHD5{1} '/data'], ...
                [channels(k) startElement(m)], [1 windowSamples]);
            end
          end
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
      obj.datatype='int16';

      
      obj = obj.getRecordingFiles(recordingFile, 'kwd');
      
      % Find the .kwe file
      %[~, name, ~] = fileparts(obj.recordingName);
      obj.fullFilename = fullfile(obj.recordingDir, obj.recordingName);
      filePrefix=strsplit(obj.recordingName,'_');
      
      triggerFile = dir([obj.recordingDir filesep filePrefix{1} '*.kwe']);
      if isempty(triggerFile)
          triggerFile = dir([obj.recordingDir filesep '*.kwe']);
          disp(['Trigger file with prefix ' filePrefix{1} ' not found!!!!, looking for other .kwe files in the same folder']);
      end
      if isempty(triggerFile)
        error('KwikRecording: Cannot file .kwe file')
      elseif length(triggerFile) > 1
        warning('KwikRecording: Multiple .kwe file found! using the first one')
      end
      obj.triggerFilename = fullfile(obj.recordingDir, triggerFile.name);
      
      if exist([obj.recordingDir filesep 'metaData.mat'],'file') && ~obj.overwriteMetaData
          obj = loadMetaData(obj); %needs recNameHD5
      else
          obj = extractMetaData(obj);
      end
      
      obj.numRecordings = length(obj.info.Groups);
      if obj.numRecordings > 1
          warning('KwikRecording: file contains multiple recordings.')
      end
      
      obj.recNameHD5 = cell(1, obj.numRecordings);
      for k = 1:obj.numRecordings
          obj.recNameHD5{k} = obj.info.Groups(k).Name;
      end
      
      obj.dataLength = obj.lengthInfo.Dataspace.Size(2);
      obj=obj.loadChLayout;
      
    end
    
    function delete(obj) %do nothing
    end
    
  end
  
  methods (Access = protected)
    
    function obj = extractMetaData(obj)
        obj.info = h5info(obj.fullFilename, obj.pathToData);
        obj.lengthInfo = h5info(obj.fullFilename, [obj.info.Groups(1).Name '/data']);
        obj.recNameHD5{1}=obj.info.Groups(1).Name;
      try
        obj.MicrovoltsPerAD = double(h5readatt(obj.fullFilename, [obj.recNameHD5{1} '/application_data'], 'channel_bit_volts'));
      catch
        obj.MicrovoltsPerAD = double(h5read(obj.fullFilename, [obj.recNameHD5{1} '/application_data/channel_bit_volts']));
      end
      
      obj.channelNumbers = 1:length(obj.MicrovoltsPerAD);
      obj.channelNames = cellfun(@(x) num2str(x), mat2cell(obj.channelNumbers,1,ones(1,numel(obj.channelNumbers))),'UniformOutput',0);
      
      try
        obj.samplingFrequency = double(h5readatt(obj.fullFilename, [obj.recNameHD5{1} '/application_data'], 'channel_sample_rates'));
      catch
        obj.samplingFrequency = double(h5read(obj.fullFilename, [obj.recNameHD5{1} '/application_data/channel_sample_rates']));
      end
      
      try
        disp('Extracting time stamp information...');
        obj.timestamps = double(h5read(obj.fullFilename, [obj.recNameHD5{1} '/application_data/timestamps']));
        disp('... done');
      catch
        disp('KwikRecording: timestamps information not available')
      end

      obj.bitDepth = double(h5readatt(obj.fullFilename, obj.recNameHD5{1}, 'bit_depth'));

      try
        obj.recordingDuration_ms = double(h5readatt(obj.fullFilename,[obj.recNameHD5{1} '/data'], 'valid_samples'));
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
      
      disp('saving meta data');
      obj.saveMetaData;
    end
  
  end
  
end

