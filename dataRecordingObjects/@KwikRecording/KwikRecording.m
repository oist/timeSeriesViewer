classdef KwikRecording < dataRecording
  %KWIKRECORDING Reads a recording in the Kwik format
  % This class is used to read kwd and kwe files containing raw data and
  % trigger information from a kwik recording. The class DOES NOT implement
  % th e whole kwik data specification but it only focuses on data and
  % triggers at the moment.
  %
  % AUTHOR: stefano.masneri@brain.mpg.de
  % DATE: 21.11.2016
  % 
  % TODO:
  %   - check how to deal with channel info for triggers
  %   - add functionalities
  
  properties
    triggerFilename; % name of the *.kwe file containing trigger info
    info;            % additional information provided in the file
    bitDepth;        % [1 x N] number of bits used to store data
    datatype;        % class of data in the recording
    sample_ms;
    fullFilename;    % path + name
  end
  
  properties (Constant = true)
    pathToData = '/recordings/';
    
    pathToTriggerData = '/event_types/TTL/events/timesamples'; % where the triggers are stored in the .kwe file
    pathToTriggerOnOff = '/event_types/TTL/events/userdata/eventID';
    pathToTriggerChannel = '/event_types/TTL/events/userdata/event_channels';
        
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
      
      windowSamples = round(window_ms / obj.sample_ms(1));
      nWindows = numel(startTime_ms);
      startTime_ms = round(startTime_ms/obj.sample_ms(1))*obj.sample_ms(1);
      startElement = round(startTime_ms/obj.sample_ms(1));
      if 0 == startElement
        startElement = 1;
      end
      %window_ms = windowSamples * obj.sample_ms(1);
      
      if isempty(channels) %if no channels are entered, get all channels
        channels=obj.channelNumbers;
      end
      nCh = numel(channels);
      
      V_uV = zeros(nCh, nWindows, windowSamples, obj.datatype); %initialize waveform matrix
      
      for k = 1:length(channels)
        V_uV(channels(k), :, :) = h5read(obj.fullFilename, [obj.channelNames{k} '/data'], ...
          [1 1 startElement], [1 nWindows windowSamples]);
      end
      
      if obj.convertData2Double
        for k = 1:size(V_uV, 1)
          V_uV(k, :, :) = double(V_uV(k, :, :)) * obj.MicrovoltsPerAD(k);
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
      
      % Assume the .kwe file has the same name as the .kwd file
      [~, name, ~] = fileparts(obj.recordingName);
      obj.fullFilename = fullfile(obj.recordingDir, obj.recordingName);
      obj.triggerFilename = fullfile(obj.recordingDir, [name, '.kwe']);
      
      fileInfo = h5info(obj.fullFilename, obj.pathToData);
      
      numCh = length(fileInfo.Groups);
      
      obj.channelNames = cell(1, numCh);
      obj.channelNumbers = zeros(1, numCh);
      obj.samplingFrequency = zeros(1, numCh);
      obj.bitDepth = zeros(1, numCh);
      obj.MicrovoltsPerAD = zeros(1, numCh);
      addMe = 0; %set to true if first channel is zero
      for k = 1:numCh
        obj.channelNames{k} = fileInfo.Groups(k).Name;
        parts = strsplit(obj.channelNames{k}, '/');
        obj.channelNumbers(k) = str2double(parts(end));
        if 1 == k
          if 0 == obj.channelNumbers(1)
            addMe = 1;
          end
        end
        obj.channelNumbers(k) = obj.channelNumbers(k) + addMe;
        obj.samplingFrequency(obj.channelNumbers(k)) = h5readatt(obj.fullFilename, ...
          obj.channelNames{k}, 'sample_rate');
        obj.bitDepth(obj.channelNumbers(k)) = h5readatt(obj.fullFilename, ...
          obj.channelNames{k}, 'bit_depth');
        tmp = h5readatt(obj.fullFilename, [obj.channelNames{k} '/application_data'], 'channel_bit_volts');
        obj.MicrovoltsPerAD(obj.channelNumbers(k)) = tmp(1);
      end
      obj.sample_ms = 1e3 / obj.samplingFrequency(1);
      
      try
        obj.recordingDuration_ms = h5readatt(obj.fullFilename, '/', 'recordingDuration');
      catch
        obj.recordingDuration_ms = 0;
      end
      
      try
        obj.datatype = h5readatt(obj.fullFilename, '/', 'datatype');
      catch
        obj.datatype = 'int16';
      end
      
    end
    
    function delete(obj) %do nothing
    end
    
  end
  
end
