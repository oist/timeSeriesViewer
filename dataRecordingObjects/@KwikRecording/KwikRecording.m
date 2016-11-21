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
  end
  
  methods
    
  end
  
  methods (Hidden = true)
    
    %class constructor
    function obj = KwikRecording(recordingFile)
      %get data files
      if nargin==0
        recordingFile=[];
      elseif nargin>1
        disp('KwikRecording: Object was not constructed since too many parameters were given at construction');
        return;
      end
      
      obj = getRecordingFile(obj, recordingFile, '.kwd');
      
    end
    
  end
  
end

