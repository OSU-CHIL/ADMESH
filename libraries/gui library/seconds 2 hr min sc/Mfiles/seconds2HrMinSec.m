function time_string = seconds2HrMinSec(time)
% seconds2HrMinSec - Convert CPU time to hours minutes seconds time string
% file.
%
% Syntax:  Open14Callback(guiFig,~,~)
%
% Inputs:
%    time - CPU time in seconds
%
% Outputs:
%    time_string - string in either hours, minutes or seconds
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dustin West
% The Ohio State University
% email address: dww.425@gmail.com
% August 2013; Last revision: 08-August-2013

%------------- BEGIN CODE --------------


if time < 60
    time_string = num2str(time,'%.2f');
    units = ' seconds.';
elseif time >= 60 && time < 3600
    time = time / 60;
    time_string = num2str(time, '%.2f');
    units = ' minutes.';
elseif time >= 3600
    time = time / 3600;
    time_string = num2str(time, '%.2f');
    units = ' hours.';
end

time_string = [time_string units];

end