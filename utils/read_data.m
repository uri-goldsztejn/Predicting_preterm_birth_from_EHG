function [X,y,name] = read_data(k,myDir)
% read_data_2 Reads the clinical information of a specified recording. 
%   C = read_data_2(k,myDir) reads sample number k in myDir.
%   Returns: X - vector of [gestational_age,age, parity, abortions, weight, hypertension, diabetes, placental_position, bleeding1, bleeding2, funneling, smoker];
%            y = vector of [is_preterm,delivery_age]; is_preterm - flag
%            indicating whether the pregnancy ended in preterm labor. 



% Check that k is valid
headerFiles = dir(fullfile(myDir,'*.hea'));
if k>length(headerFiles)
    error('k exceeds number of files\n');
end


% Get filename
signal_baseFileName = headerFiles(k).name;
name = signal_baseFileName;
if isequal(signal_baseFileName(end-9),'_')
    header_fullFileName = fullfile(myDir, [signal_baseFileName(1:end-5),'.hea']);
else
    header_fullFileName = fullfile(myDir, [signal_baseFileName(1:end-4),'.hea']);
end


%get days to labor
[is_preterm,gestational_age,delivery_age] = get_preterm_flag(header_fullFileName);

% Get age
buffer = fileread(header_fullFileName) ;
substr = 'Age';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    age = nan;
else
    age_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    
    if strcmp(age_str,"no") || strcmp(age_str,"No") || strcmp(age_str,"none") || strcmp(age_str,"None")...
            || strcmp(age_str,"negative") || strcmp(age_str,"Negative")
        age = nan;
    else
        age = str2num(age_str);
    end
end

% Get parity
buffer = fileread(header_fullFileName) ;
substr = 'Parity';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    parity =nan;
else
    
    parity_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    if strcmp(parity_str,"no") || strcmp(parity_str,"No") || strcmp(parity_str,"none") || strcmp(parity_str,"None")...
            || strcmp(parity_str,"negative") || strcmp(parity_str,"Negative")
        parity = 0;
    else
        parity = str2num(parity_str);
    end
end

% Get number of past abortions
buffer = fileread(header_fullFileName) ;
substr = 'Abortions';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    abortions =nan;
else
    abortions_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+7), '%s', 1);
    
    if strcmp(abortions_str,"no") || strcmp(abortions_str,"No") || strcmp(abortions_str,"none") || strcmp(abortions_str,"None")...
            || strcmp(abortions_str,"negative") || strcmp(abortions_str,"Negative")
        abortions = 0;
    else
        abortions = str2num(abortions_str);
    end
end

% Get weight
buffer = fileread(header_fullFileName) ;
substr = 'Weight';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    weight =nan;
else
    weight_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    if strcmp(weight_str,"no") || strcmp(weight_str,"No") || strcmp(weight_str,"none") || strcmp(weight_str,"None")...
            || strcmp(weight_str,"negative") || strcmp(weight_str,"Negative")
        weight = nan;
    else
        weight = str2num(weight_str);
    end
    
end

% Get hypertension
buffer = fileread(header_fullFileName) ;
substr = 'Hypertension';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    hypertension =nan;
else
    hypertension_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    if strcmp(hypertension_str,"no") || strcmp(hypertension_str,"No") || strcmp(hypertension_str,"none") || strcmp(hypertension_str,"None")...
            || strcmp(hypertension_str,"negative") || strcmp(hypertension_str,"Negative")
        
        hypertension = 0;
    else
        hypertension = 1;
    end
    
end


% Get diabetes
buffer = fileread(header_fullFileName) ;
substr = 'Diabetes';
loc    = strfind(buffer, substr) ;

if isempty(loc)
    diabetes =nan;
else  
    diabetes_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    if strcmp(diabetes_str,"no") || strcmp(diabetes_str,"No") || strcmp(diabetes_str,"none") || strcmp(diabetes_str,"None")...
            || strcmp(diabetes_str,"negative") || strcmp(diabetes_str,"Negative")
        diabetes = 0;
    else
        diabetes = 1;
    end
    
end

% Get Placental_position
buffer = fileread(header_fullFileName) ;
substr = 'Placental_position';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    placental_position =nan;
else
    placental_position_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    
    if strcmp(placental_position_str,"no") || strcmp(placental_position_str,"No") || strcmp(placental_position_str,"none") || strcmp(placental_position_str,"None")...
            || strcmp(placental_position_str,"negative") || strcmp(placental_position_str,"Negative")
        
        placental_position = nan;     
    elseif strcmp(placental_position_str,"front")
        placental_position = 1;
    elseif strcmp(placental_position_str,"end")
        placental_position = 0;       
    end
end


% Get Bleeding_first_trimester
buffer = fileread(header_fullFileName) ;
substr = 'Bleeding_first_trimester';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    bleeding1 = nan;
else   
    bleeding1_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    if strcmp(bleeding1_str,"no") || strcmp(bleeding1_str,"No") || strcmp(bleeding1_str,"none") || strcmp(bleeding1_str,"None")...
            || strcmp(bleeding1_str,"negative") || strcmp(bleeding1_str,"Negative")  
        bleeding1 = 0;
    else
        bleeding1 = 1;
    end
    
end

% Get Bleeding_second_trimester
buffer = fileread(header_fullFileName) ;
substr = 'Bleeding_second_trimester';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    bleeding2 =nan;
else
    bleeding2_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    if strcmp(bleeding2_str,"no") || strcmp(bleeding2_str,"No") || strcmp(bleeding2_str,"none") || strcmp(bleeding2_str,"None")...
            || strcmp(bleeding2_str,"negative") || strcmp(bleeding2_str,"Negative") 
        bleeding2 = 0;
    else
        bleeding2 = 1;
    end
    
end

% Get funneling
buffer = fileread(header_fullFileName) ;
substr = 'Funneling';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    funneling =nan;
else
    funneling_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+8), '%s', 1);
    if strcmp(funneling_str,"no") || strcmp(funneling_str,"No") || strcmp(funneling_str,"none") || strcmp(funneling_str,"None")...
            || strcmp(funneling_str,"negative") || strcmp(funneling_str,"Negative")
        
        funneling = 0;
    else
        funneling = 1;
    end
    
end


% Get smoker
buffer = fileread(header_fullFileName) ;
substr = 'Smoker';
loc    = strfind(buffer, substr) ;
if isempty(loc)
    smoker =nan;
else
    smoker_str   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+2), '%s', 1);
    if strcmp(smoker_str,"no") || strcmp(smoker_str,"No") || strcmp(smoker_str,"none") || strcmp(smoker_str,"None")...
            || strcmp(smoker_str,"ne") || strcmp(smoker_str,"Ne")        
        smoker = 0;
    else
        smoker = 1;
    end
    
end

y = [is_preterm,delivery_age];
X = [gestational_age,age, parity, abortions, weight,...
    hypertension, diabetes, placental_position, bleeding1, bleeding2, funneling, smoker];
end