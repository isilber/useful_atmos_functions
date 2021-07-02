%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Israel Silber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sonde = Get_sonde_data(Station_id, 00Z_or_12Z_flag, wanted_datenum)
% Sonde = Get_sonde_data(Station_id, 00Z_or_12Z_flag, wanted_year, wanted_month, wanted_day)
% The function downloads radiosonde data from the Wyoming website, and
% saves it in the 'Sonde' structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs include:
% 'Station_id' - the number id of the wanted station, e.g., Green Bay
% has id #72645.
% '00Z_or_12Z_flag' - determine if 00Z or 12Z should be acquired - true for
% 12Z and false for 00Z.
% The wanted day should be sent to the function via a single datenum
% variable or a (wanted_year, wanted_month, wanted_day) double sequence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: all inputs should be of type double.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function Sonde = Get_sonde_data(Station_id, Z12_true_Z00_false, varargin)
function Sonde = Get_sonde_data(Station_id, Hour, varargin)

%% used for debugging.
% Station_id = 72645; % 72645 is for Wisconsin.
% Year = '2016';
% Month = '01';
% Day = '28';
% Z00_true_Z12_false = true;  % true = extracting 00Z, false = extracting 12Z;

%% Determine wanted day and time from input.
if nargin == 2 || nargin == 4
    disp('ERROR: not enough input arguments')
    Sonde = nan;
    return
elseif nargin == 3
    Year = datestr(varargin{1}, 'yyyy');
    Month = datestr(varargin{1}, 'mm');
    Day = datestr(varargin{1}, 'dd');
elseif nargin > 3
    Year = num2str(varargin{1});
    Month = num2str(varargin{2});
    Day = num2str(varargin{3});
end

% if Z12_true_Z00_false
%     Hour = '12';
% else
%     Hour = '00';
% end

%% retrieve data.
file_str = urlread(['http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=', Year, '&MONTH=', Month, '&FROM=', Day, Hour, '&TO=', Day, Hour, '&STNM=',num2str(Station_id)]);
valid_data = strfind(file_str, ['Can''t get ', num2str(Station_id)]); % searching for an error message (missing data) in the retrieved sonde data. If empty, then data exists...

% arranging data only if avalilable.
if isempty(valid_data)
    
    start_ind = strfind(file_str, '-----------------------------------------------------------------------------'); % this string is 77 characters in length.
    end_ind = strfind(file_str, 'Station information');
    file_str = file_str(start_ind(end) + 78: end_ind); % 77 + 1 to include the '\n' character and start from the first pressure field.
    counter = 1; % count number of lines.
    sonde_tmp = NaN(floor(floor(length(file_str)/7)/ 11), 11); % radiosonde has 11 fields each has a length of 7 characters.
    for ii = 1: floor(length(file_str)/7)
        if ii * 7 + counter < length(file_str)
            strtmp = file_str((ii-1) * 7 + 1 + counter: ii * 7 + counter); % grabbing a single parameter a time.
            find_digit = find(isstrprop(strtmp, 'digit') | isstrprop(strtmp, 'punct') );
            if ~isempty(find_digit) % some recorded parameters are sometimes missing.
                sonde_tmp(counter, mod(ii, 11) + (mod(ii, 11) == 0)*11) = str2double(strtmp(find_digit(1): find_digit(end)));
            end
            if mod(ii, 11) == 0 % 11 variable, so looking for the last one.
                counter = counter + 1;
            end
        end
    end
    sonde_tmp = sonde_tmp(1: end-1, :); % delete the last line which contains only nans.
    
    %% arrange data in a single struct.
    Sonde.pressure = sonde_tmp(:, 1); % pressure in mB
    Sonde.alt = sonde_tmp(:, 2); % altitude in m.
    Sonde.drybulb_temp = sonde_tmp(:, 3); % Air temperature in C.
    Sonde.dewpoint_temp = sonde_tmp(:, 4); % Dew point temperature in C.
    Sonde.RH = sonde_tmp(:, 5);% RH in %.
    Sonde.mix_ratio = sonde_tmp(:, 6); % Mixing ratio in g/kg.
    Sonde.wind_direction = sonde_tmp(:, 7); % in direction in degrees clockwise from north.
    Sonde.wind_speed = sonde_tmp(:, 8).* 0.514444444; % wind speed in m/s (after multiplying with the conversion factor).
    Sonde.theta = sonde_tmp(:, 9); % potential temperature in K.
    Sonde.theta_e = sonde_tmp(:, 10); % equivalent potential temperature in K.
    % Sonde.theta_tv = sonde_tmp(:, 11); % virtual potential temperature in K.
    
else % if no data exists, Sonde is returned empty.
    Sonde = [];
    disp(['missing data for ', Year, Month, Day])
end