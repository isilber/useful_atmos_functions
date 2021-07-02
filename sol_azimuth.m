%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Israel Silber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculationg solar azimuthn angle. May return a vector or scalar -
%according to the 'day' parameter.
%aveH is the UT hour on  a double format.- thus 12:30 will become 12.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [solAZ]=sol_azimuth(lat,lon,Day)

% datevc of the wanted time.
Dayvec = datevec(Day);
aveH = Dayvec(:, 4) + Dayvec(:, 5)./ 60 + Dayvec(:, 6)./ 3600;
Day = datenum(Dayvec(:,1), Dayvec(:,2), Dayvec(:,3));
clear Dayvec

%calculating days since January  1st of the year.
N=Day-datenum(str2double(datestr(Day(1),10)),1,1);

%solar local hour angle, aveH is the UT hour
solH=15*aveH-180+lon;
solH(solH > 180) = solH(solH > 180) - 360; % correcting for longitude issues.
%solar declenation - accurate using Earth's orbit parameters.
%see: http://en.wikipedia.org/wiki/Position_of_the_Sun
solD=asind(sind(-23.44).*cosd(360./365.24.*(N+10)+360./pi.*0.0167.*sind(360./365.24.*(N-2))));

%solar elevation angle calculation.
%see: http://en.wikipedia.org/wiki/Solar_elevation_angle
SZA=90 - asind(cosd(solH).*cosd(solD).*cosd(lat)+sind(solD).*sind(lat));

solAZ = acosd((sind(solD).*cosd(lat) - cosd(solH).*cosd(solD).*sind(lat))./ sind(SZA));
solAZ(solH > 0) = 360 - solAZ(solH > 0); % correction due to the cos symmetry around 0.
% solAZ = asind(-sind(solH).* cosd(solD)./ sind(SZA)); % alternative

return

