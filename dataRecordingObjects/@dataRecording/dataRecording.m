classdef (Abstract) dataRecording < handle
    properties
        recordingName %(String) The name of the recording
        recordingDir % (String) Full directory containing the recorded session
        dataFileNames % (Cell 1 x N)  array of N recording data file names
        startDate %(1x1) Start date (time) of Recording (matlab long format)
        endDate %(1x1) End date (time) of Recording (matlab long format)
        samplingFrequency %(1x1) Sampling rate
        recordingDuration_ms %(1x1) the total duration of the recording in [ms]
        channelNames % (Cell 1xN) a cell array with the N names of the channels
        channelNumbers % (1xN) an array with integer channel numbers
        channelNumbersOrignal
        triggerNames %the names of trigger channels
        dspLowCutFrequency % (1x1) Low-pass cutoff frequency in the Neuralynx DSP (in raw data)
        dspHighCutFrequency % (1x1) High-pass cutoff frequency in the Neuralynx DSP (in raw data)
        nRecordings % (1x1) number of recording files
        chLayoutNumbers %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid
        chLayoutNames %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid
        chLayoutPositions % (1xN or 2xN or 3xN) array of electrode position in [x or x,y or x,y,z]  
        layoutName %the name of the channel layout (electrode type)
        n2s % a translation between the number of the channel to the serial number of the channel (in the case where all channels are consecutive)
        
        convertData2Double = 1 % if data should be converted to double from the original quantization
        ZeroADValue
        MicrovoltsPerAD
        overwriteMetaData = false;
    end
    
    properties (SetAccess=protected) %these are properties that are not synchronized or loaded from meta files
        multifileMode %(logical 1x1) if multi files were selected
        folderMode
    end
    
    properties (Constant, Abstract)
        defaultLocalDir %Default directory from which search starts
        signalBits % the quantization of the sampling card
        numberOfCharFromEndToBaseName %the number of characters from the end of the file name to get to the base file name
    end
    methods
        function delete(obj) %closing all open files when object is deleted
            obj=closeOpenFiles(obj);
        end
        function obj=closeOpenFiles(obj)
        end
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the recording (if empty takes the default name)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
        end
        function [V_uV,T_ms]=getAnalogData(obj,channels,startTime_ms,window_ms,name)
            %Extract recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getAnalogData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the recording (if empty takes the default name)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
        end
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,endTime_ms,TTLbits)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name - the name of the trigger (if empty takes the default name)
            %Output: T_ms - trigger times [ms]
        end
        function [D,T_ms]=getDigitalData(obj,startTime_ms,window_ms,name)
            %Extract MCRack digital data from file to memory
            %Usage: [D,T_ms]=getDigitalData(startTime_ms,window_ms,name)
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: D - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
        end
        function saveMetaData(obj)
            %Save object properties (metaData) to file
            %Usage : obj.saveMetaData;
            %Input : []
            props.metaClassData=metaclass(obj);
            props.allPropName={props.metaClassData.PropertyList.Name}';
            props.allPropIsConstant=cell2mat({props.metaClassData.PropertyList.Constant}');
            props.allPropSetAccess={props.metaClassData.PropertyList.SetAccess}';

            pNonConstantProps=find(~props.allPropIsConstant & ~strcmp(props.allPropSetAccess,'protected'));
            for i=1:numel(pNonConstantProps)
                metaData.(props.allPropName{pNonConstantProps(i)})=obj.(props.allPropName{pNonConstantProps(i)});
            end
            save([obj.recordingDir filesep 'metaData'],'metaData');
        end
        
        function deleteMetaData(obj)
            if ~iscell(obj.recordingDir)
                delete([obj.recordingDir filesep 'metaData.mat']);
            else
                for i=1:numel(obj.recordingDir)
                    delete([obj.recordingDir{i} filesep 'metaData.mat']);
                end
            end
        end
        
        function obj=loadMetaData(obj,fileName)
            %Load object properties (metaData) from file
            %Usage : obj.loadMetaData;
            %Input : fileName - if entered loads meta data from this file, else loads data from main data directory
            oldRecordingDir=obj.recordingDir;
            if ~iscell(obj.recordingDir) %regular recording
                if nargin==2
                    load(fileName);
                else
                    load([obj.recordingDir filesep 'metaData'],'metaData');
                end
                fieldNames=fieldnames(metaData);
                for i=1:numel(fieldNames)
                    obj.(fieldNames{i})=metaData.(fieldNames{i});
                end
            else %multi file recording 
                for i=1:numel(obj.recordingDir)
                    if nargin==2
                        load(fileName{i});
                    else
                        load([oldRecordingDir{i} filesep 'metaData'],'metaData');
                    end
                    fieldNames=fieldnames(metaData);
                    for j=1:numel(fieldNames)
                        if numel(metaData.(fieldNames{j}))==1
                            obj.(fieldNames{j})(i)=metaData.(fieldNames{j});
                        else
                            if ~iscell(obj.(fieldNames{j}))==1
                                obj.(fieldNames{j})=cell(1,numel(obj.recordingDir));
                            end
                            obj.(fieldNames{j}){i}=metaData.(fieldNames{j});
                        end
                    end
                end
            end
            obj.recordingDir=oldRecordingDir;
        end
        function obj=loadChLayout(obj)
            %checks for a .chMap file with the recording name in the same folder of the recording and extract the layout information
            %txt should correspond to layout file name on path
            if iscell(obj.recordingDir)
                recordingDir=obj.recordingDir{1};
            end
            chMapFiles=dir([obj.recordingDir filesep '*.chMap']);
            chMapFiles={chMapFiles.name};
            switch numel(chMapFiles)
                case 0 %not channel map file found
                    disp('No .chMap files were found for this recording');
                    return;
                case 1 %there is only one channel map file, this file will apply to all the recordings
                    chMapFiles=[obj.recordingDir filesep chMapFiles{1}];
                otherwise %there are several files, in which case each recording should have its own channel map file with the appropriate name
                    chMapFiles=dir([obj.recordingDir filesep obj.recordingName(1:end-numel(obj.fileExtension)-1) '*.chMap']);
                    chMapFiles={chMapFiles.name};
                    if numel(chMapFiles)~=1
                        disp('Channel map file name (*.chMap) does not correspond to the recording file name');
                        return;
                    else
                        chMapFiles=[obj.recordingDir filesep chMapFiles{1}];
                    end
            end
            
            disp(['channel map extracted from ' chMapFiles]);
            A = importdata(chMapFiles);
            if numel(A)==1
                obj.layoutName=['layout_' A{1}];
                load(obj.layoutName);
                obj.chLayoutNumbers=En;
                obj.chLayoutNames=Ena;
                obj.chLayoutPositions=Enp;
            else
                for i=1:numel(A)
                    obj.layoutName{i}=['layout_' A{i}];
                    load(obj.layoutName{i});
                    obj.chLayoutNumbers{i}=En;
                    obj.chLayoutNames{i}=Ena;
                    obj.chLayoutPositions{i}=Enp;
                end
                
            end
            
        end
        
        function obj=getRecordingFiles(obj,recordingFile,fileExtension)
            %Get directory with data files 
            %Usage: obj = getRecordingFiles(obj,recordingFile,fileExtension)
            %if no recording file is entered lauches a GUI
            %if no file extension entered, a directory is chosen rather than a specific files (for example for neuralynx recordings)
            
            %If no files were entered open GUI for choosing a file or a directory else get the files entered
            if ~isempty(recordingFile) %if directory with data was not entered open get directory GUI
                
                obj.multifileMode=iscell(recordingFile);
                if obj.multifileMode,singleRecordingFile=recordingFile{1};else singleRecordingFile=recordingFile;end
                if isdir(singleRecordingFile)
                    obj.folderMode=true; %a folder is chosen and the files inside examined
                else
                    obj.folderMode=false; %a file or list of files is selected
                end
                
                if ~obj.folderMode
                    if ~obj.multifileMode
                        recordingFile={recordingFile};
                    end
                    obj.nRecordings=numel(recordingFile);
                    for i=1:obj.nRecordings
                        [pathstr, name{i}, ext] = fileparts(recordingFile{i});
                        obj.dataFileNames{i}=[name{i} ext];
                        if ~exist([pathstr filesep obj.dataFileNames{i}],'file')
                            disp(['Searching for recording file: ' [pathstr filesep obj.dataFileNames{i}]]);
                            error('Object was not constructed since no valid file name was chosen');
                        end
                    end
                else
                    if ~obj.multifileMode
                        [pathstr, name] = fileparts(recordingFile);
                        obj.dataFileNames{1}=recordingFile;
                    else
                        [pathstr, name] = cellfun(@(x) fileparts(x),recordingFile,'UniformOutput',0);
                        obj.dataFileNames=recordingFile;
                    end
                end
                
                if isempty(pathstr) %in case the file is in the current directory
                    if ispc
                        obj.recordingDir=[cd filesep];
                    end
                else
                    obj.recordingDir=pathstr;
                    if ispc
                        if ~iscell(obj.recordingDir)
                            obj.recordingDir=[obj.recordingDir filesep];
                        else
                            obj.recordingDir=cellfun(@(x) [x filesep],obj.recordingDir,'UniformOutput',0);
                        end
                    end
                end
                
                if obj.multifileMode & obj.folderMode %some of the condition below can be removed
                    if ~isdir(obj.recordingDir{1})
                        error('Object was not constructed since no valid folder was choosen');
                    end
                else
                    if ~isdir(obj.recordingDir)
                        error('Object was not constructed since no valid folder was choosen');
                    end
                end
            else %if directory with data was not entered open get directory GUI
                if ~obj.folderMode
                    [obj.dataFileNames,obj.recordingDir]= uigetfile(['*.' fileExtension],['Choose the ' fileExtension ' file'],obj.defaultLocalDir,'MultiSelect','on');
                    if ~iscell(obj.dataFileNames)
                        obj.dataFileNames={obj.dataFileNames};
                    end
                    if obj.dataFileNames{1}==0 %no folder chosen
                        disp('Object was not constructed since no folder was choosen');
                        return;
                    end
                    obj.nRecordings=numel(obj.dataFileNames);
                    if obj.nRecordings>1
                        obj.multifileMode=true;
                    end
                else
                    [obj.recordingDir]= uigetdir(obj.defaultLocalDir,'Choose the data folder');
                    [pathstr, name] = fileparts(obj.recordingDir);
                    obj.recordingDir=[pathstr filesep];
                    obj.multifileMode=false;
                end
            end
            if ~obj.folderMode
                obj.recordingName=obj.dataFileNames{1};
            else
                obj.recordingName=name;
            end
        end
    end
end