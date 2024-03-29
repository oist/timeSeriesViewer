classdef (Abstract) dataRecording < handle
    properties
        recordingName %(String) The name of the recording
        recordingDir % (String) Full directory containing the recorded session
        dataFileNames % (Cell 1 x N)  array of N recording data file names
        startDate %(1x1) Start date (time) of Recording (matlab long format)
        endDate %(1x1) End date (time) of Recording (matlab long format)
        samplingFrequency %(1xN) Sampling rate
        recordingDuration_ms %(1x1) the total duration of the recording in [ms]
        channelNames % (Cell 1xN) a cell array with the N names of the channels
        channelNumbers % (1xN) an array with integer channel numbers
        channelNumbersOrignal
        triggerNames %the names of trigger channels
        analogChannelNumbers
        dspLowCutFrequency % (1x1) Low-pass cutoff frequency in the Neuralynx DSP (in raw data)
        dspHighCutFrequency % (1x1) High-pass cutoff frequency in the Neuralynx DSP (in raw data)
        nRecordings % (1x1) number of recording files
        chLayoutNumbers %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid
        chLayoutNames %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid
        chLayoutPositions % (1xN or 2xN or 3xN) array of electrode position in [x or x,y or x,y,z]
        layoutName %the name of the channel layout (electrode type)
        n2s % a translation between the number of the channel to the serial number of the channel (in the case where all channels are consecutive)
        
        convertData2Double = 1; % if data should be converted to double from the original quantization
        ZeroADValue
        MicrovoltsPerAD
        datatype        % class of data in the recording
        
        overwriteMetaData = false;
        electrodePitch
        metaDataFile
    end
    
    properties (SetAccess=protected) %these are properties that are not synchronized or loaded from meta files
        multifileMode %(logical 1x1) if multi files were selected
        folderMode = false;
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
%            save([obj.recordingDir filesep obj.recordingName '_metaData.mat'],'metaData');
        end
        
        function [X,Y,Z]=getElectrodePositions(obj,electrodePitch)
            %if recording object contains electrode positions, use these, if not
            if ~isempty(obj.chLayoutPositions)
                disp('Getting positions from layout files');
                X=obj.chLayoutPositions(1,:);
                Y=obj.chLayoutPositions(2,:);
            else
                if nargin==2
                    obj.electrodePitch=electrodePitch;
                elseif nargin==1 & isempty(obj.electrodePitch)
                    obj.electrodePitch=100;
                end
                disp(['Getting positions from grid layout, assuming pitch of ' num2str(obj.electrodePitch) 'um !!!!!']);
                
                %Build inverse map between electrode and location
                [meshX,meshY]=meshgrid(1:size(obj.chLayoutNumbers,1),1:size(obj.chLayoutNumbers,2));
                X(obj.chLayoutNumbers(~isnan(obj.chLayoutNumbers)))=meshX(~isnan(obj.chLayoutNumbers))*obj.electrodePitch;
                Y(obj.chLayoutNumbers(~isnan(obj.chLayoutNumbers)))=meshY(~isnan(obj.chLayoutNumbers))*obj.electrodePitch;
            end
            Z=zeros(size(Y));
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
            try
                oldRecordingDir=obj.recordingDir;
                if ~iscell(obj.recordingDir) %regular recording
                    if nargin==2
                        load(fileName);
                    else
                        load([obj.recordingDir filesep obj.recordingName '_metaData'],'metaData');
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
                            load([oldRecordingDir{i} filesep obj.recordingName '_metaData'],'metaData');
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
            catch errorMsg
                disp('Error while extracting fields from meta data. Trying re-extract meta data...');
                obj=extractMetaData(obj);
            end
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
            
            A = importdata(chMapFiles);
            elecString=regexp(A{1},'_','split');
            obj.electrodePitch=str2num(elecString{1});
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
            fprintf('Channel map with pitch %d and layout %s extracted from %s\n',obj.electrodePitch,elecString{2},chMapFiles);
            
            %check that all recorded channels are contained within the layout
            if numel(obj.channelNumbers)>numel(intersect(obj.channelNumbers,En(:)))
                warning('Notice that some of the recorded channels are not contained in the layout file, this may result in errors in further analysis!');
            end
            
        end
        
        function []=convertLayoutKSort(obj,outputFile,badChannels)
            if nargin<2
                outputFile=fullfile(obj.recordingDir, 'chanMap.mat');
            end
            if nargin<3
                badChannels=[];
            end
            
            % here I know a priori what order my channels are in.  So I just manually
            % make a list of channel indices (and give an index to dead channels too). chanMap(1) is the row in the raw binary
            % file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to be dead channels.
            
            chanMap = obj.channelNumbers;
            
            % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
            % Now we declare which channels are "connected" in this normal ordering,
            % meaning not dead or used for non-ephys data
            
            badChannels=[];
            connected = true(numel(chanMap), 1);
            if ~isempty(badChannels)
                connected(badChannels) = false;
            end
            
            % now we define the horizontal (x) and vertical (y) coordinates of these
            % 34 channels. For dead or nonephys channels the values won't matter. Again
            % I will take this information from the specifications of the probe. These
            % are in um here, but the absolute scaling doesn't really matter in the
            % algorithm.
            
            xcoords = obj.chLayoutPositions(1,:);
            ycoords = obj.chLayoutPositions(2,:);
            
            % Often, multi-shank probes or tetrodes will be organized into groups of
            % channels that cannot possibly share spikes with the rest of the probe. This helps
            % the algorithm discard noisy templates shared across groups. In
            % this case, we set kcoords to indicate which group the channel belongs to.
            % In our case all channels are on the same shank in a single group so we
            % assign them all to group 1.
            
            kcoords = true(numel(chanMap), 1);
            
            % at this point in Kilosort we do data = data(connected, :), ycoords =
            % ycoords(connected), xcoords = xcoords(connected) and kcoords =
            % kcoords(connected) and no more channel map information is needed (in particular
            % no "adjacency graphs" like in KlustaKwik).
            % Now we can save our channel map for the eMouse.
            
            % would be good to also save the sampling frequency here
            fs = obj.samplingFrequency;
            
            save(outputFile, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
            
            disp(['Channel information saved in : ' outputFile]);
            
        end

        function []=convertLayouteJRClust(obj,padSize,outputName)
            %convertLayouteJRClust(obj,padSize,outputName)
            %Make probe (.prb) file for using with jrclust
            %pad size - [height (y),widht (x)]
            %outputName - a name of the prb file (e.g. probe1.prb)        nCh=numel(~isnan(Enp(1,:)));
            %fid=fopen('layout_100_12x12.prb','w');
            if nargin<3
                outputName=[obj.recordingDir filesep obj.layoutName '_JRC.prb'];
            end
            if nargin<2
                error('Pad size must be entered as an external parameter');
            end
            fid=fopen(outputName,'w');
            
            nCh=size(obj.chLayoutPositions,2);
            fprintf(fid, 'channels = [1:%d];\n\n',nCh);
            fprintf(fid, 'geometry = [%.1f,%.1f',obj.chLayoutPositions(1,1),obj.chLayoutPositions(2,1));
            for i=2:nCh
                fprintf(fid,';%.1f,%.1f',obj.chLayoutPositions(1,i),obj.chLayoutPositions(2,i));
            end
            fprintf(fid, '];\n\n');
            
            fprintf(fid, 'pad = [%.1f,%.1f];\n\n',padSize(1),padSize(2));
            
            fprintf(fid, 'cviShank = {1:%d};',nCh);
            fclose(fid);
        end
        
        function export2Binary(obj,targetFile,dataChannels,medianFilterGroup)
            tic;
            targetDataType='int16';
            if nargin<4
                medianFilterGroup=[];
            else
                for i=1:numel(medianFilterGroup)
                    [~,pGroup{i}]=intersect(dataChannels,medianFilterGroup{i});
                end
            end
            
            if ~strcmp(obj.datatype,targetDataType)
                fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!',obj.datatype,targetDataType);
                zeroValue=2^16/2;
                convertDataType=true;
            else
                zeroValue=0;
                convertDataType=false;
            end
            
            %converts data recording object to kilo sort binary format for sorting
            if nargin<3
                dataChannels=obj.channelNumbers;
            end
            if nargin<2 || isempty(targetFile)
                targetFile=[obj.recordingDir filesep 'binaryData.bin'];
            end
            if ~any(strcmp(targetFile(end-3:end),{'.dat','.bin'}))
                error('input file should have a ''.dat/.bin'' extension');
            end
            
            chunkSize=2*60*1000; %msec
            startTimes=0:chunkSize:obj.recordingDuration_ms;
            endTimes=[startTimes(2:end) obj.recordingDuration_ms];
            
            if ~exist(targetFile,'file')
                %open data file
                fid = fopen(targetFile, 'w+');
                obj.convertData2Double=0;
                
                fprintf('\nConverting blocks to binary %s format(/%d) : ',targetDataType,numel(startTimes));
                nDigits=0;
                for j=1:numel(startTimes)
                    fprintf([repmat('\b',1,nDigits) '%d'],j);nDigits=length(num2str(j));
                    if ~convertDataType
                        data=squeeze(obj.getData(dataChannels,startTimes(j),endTimes(j)-startTimes(j)));
                    else
                        data=int16(int32(squeeze(obj.getData(dataChannels,startTimes(j),endTimes(j)-startTimes(j))))-zeroValue);
                    end
                    if ~isempty(medianFilterGroup)
                        for i=1:numel(medianFilterGroup)
                            data(pGroup{i},:)=bsxfun(@minus,data(pGroup{i},:),median(data(pGroup{i},:)));
                        end
                    end
                    pause(0.001);
                    fwrite(fid, data, ['*' targetDataType]);
                end
                fclose(fid);
                fprintf('\nConversion complete (binary %s)\n',targetDataType);
            else
                disp('file already exists, please delete data first and run again!');
            end
            
            fprintf('Writing trigger file...\n');
            T=obj.getTrigger;
            nT=cellfun(@(x) numel(x),T);
            pT=find(nT>0);
            
            triggerFile=[targetFile(1:end-4) '_Triggers.bin'];
            fid = fopen(triggerFile, 'w+');
            fwrite(fid,uint32(nT+1),'*uint32');
            for i=1:numel(pT)
                fwrite(fid, uint32(T{pT(i)}*obj.samplingFrequency(1)/1000)+1,'*uint32');
            end
            fclose(fid);
            
            metaDataFile=[targetFile(1:end-4) '.meta'];
            if ~exist(metaDataFile,'file')
                fid=fopen(metaDataFile,'w');
                fprintf(fid,'nSavedChans = %d\n',numel(dataChannels));
                fprintf(fid,'sRateHz = %d\n',obj.samplingFrequency(1));
                fprintf(fid,'nChans = %d\n',numel(dataChannels));
                fprintf(fid,'nTriggerChans = %d\n',numel(nT));
                fprintf(fid,'nAnalogChans = %d\n',numel(obj.analogChannelNumbers));
                fprintf(fid,'vcDataType = %s\n',targetDataType);
                fprintf(fid,'scale = %.12f\n',obj.MicrovoltsPerAD(1));
                fprintf(fid,'vcProbe = %s\n',obj.layoutName);
                fclose(fid);
            else
                disp(['Meta data file: ' metaDataFile ' alreday exists!, please first delete']);
            end
            toc;
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
                            error('Object was not constructed since no valid recording file name was chosen');
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
                [~,obj.recordingName]=fileparts(obj.dataFileNames{1});
            else
                obj.recordingName=name;
            end
            obj.metaDataFile=[obj.recordingDir filesep obj.recordingName '_metaData'];
        end
    end
end