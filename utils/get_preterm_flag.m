function [is_preterm,gestational_age,delivery_age] = get_preterm_flag(header_fullFileName)
% Reads a header file and returns a flag indicating whther the labor was premature, the gestational age at recording, 
% and the remaining days to labor


buffer = fileread(header_fullFileName) ;
substr = 'Rectime';
loc    = strfind(buffer, substr) ;
gestational_age_string   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+4), '%s', 1);
if length(gestational_age_string)<3
    gestational_age = str2num(gestational_age_string(1:2))*7 ;
else
    gestational_age = str2num(gestational_age_string(1:2))*7 + round(str2num(gestational_age_string(4))/10*7);
end

substr = 'Gestation';
loc    = strfind(buffer, substr) ;
delivery_age_string   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+4), '%s', 1);

if(isempty(delivery_age_string))
    substr = 'age at delivery(w/d):';
    loc    = strfind(buffer, substr) ;
    delivery_age_string   = sscanf(buffer(loc+numel(substr):loc+numel(substr)+4), '%s', 1);
end

if length(delivery_age_string)<3
    delivery_age = str2num(delivery_age_string(1:2))*7 ;
else
    delivery_age = str2num(delivery_age_string(1:2))*7 + round(str2num(delivery_age_string(4))/10*7);
end

% preterm flag
if delivery_age < 259
    is_preterm = 1;
else
    is_preterm = 0;
end


end

