%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Israel Silber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Output_struct] = calculate_theta_and_more(T, p, RH)
% [Output_struct] = calculate_theta_and_more(T, p, '', w)
% [Output_struct] = calculate_theta_and_more(__, model)
% [Output_struct] = calculate_theta_and_more(__, model, Uncertainties)
% The function receives three parameters and calculates others.
% Useful for working with radiosonde data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% T - temperatures in Celsius (can be an array).
% p - pressure in mbar (can be an array).
% RH - relative humidity in % (can be an array).
% water_w - water vapor mixing ratio [g/kg] (input 4 - if exists, ignoring RH).
% Uncertainties - two component array with T and RH uncertainties,
% respectively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output stucture contains:
% L_v - latent heat of vaporization (as it is temperature dependent).
% L_s - latent heat of sublimation (as it is temperature dependent).
% e_s - water vapor saturation pressure above liquid water.
% e_i - water vapor saturation pressure above ice water.
% e - water vapor pressure.
% w - water vapor mixing ratio.
% Theta - potential temperature.
% Theta_e - equivalent potential temperature.
% Theta_v - virtual potential temperature.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to define the model without w (with RH), use an empty variable instead of
% w.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Output_struct] = calculate_theta_and_more(T, p, RH, varargin)

if nargin >= 4 % needed for the retrieval of RH from w.
    if ~isempty(varargin{1}) % make sure not to use w accidentally (if only the model was needed to be specified).
        water_w = varargin{1}./1000;
    end
    if nargin >= 6 % needed for theta_e uncertainty calculation.
        if ~isempty(varargin{3})
            delta_Tk = varargin{3}(1);
            delta_RH = varargin{3}(2);
        end
    end
end
p = p.*100;
%% Determine constants and reference values.
p_0 = 1e5; % mBar or hPa (using only as ratio so mb is ok).
R_d = 287.058; % gas constant for dry air [J/(kg*K)].
R_v = 461.5; % gas constant for water vapor [J/(kg*K)].
c_p = 1005.7; % +- 2.5 [J/(kg*K)] - specific heat capacity of dry air at 273K in a constant pressure.

% liquid water saturations polynom coefficients based on Flatau et al., 1992.
p_water_s_coeff = [0.209339997e-15...
    -0.373208410e-12...
    0.892344772e-10...
    0.196237241e-7...
    0.305903558e-5...
    0.264461437e-3...
    0.143064234e-1...
    0.444007856...
    6.11213476];

% ice water saturations polynom coefficients based on Flatau et al., 1992.
p_water_i_coeff = [0.262655803e-14...
    0.149436277e-11...
    0.387940929e-9...
    0.602780717e-7...
    0.614396778e-5...
    0.420547422e-3...
    0.188369801e-1...
    0.503109514...
    6.11123516];

%% Control variables.
if nargin >= 5 % specify model
    if ~isempty(varargin{2})
        use_Flatau = varargin{2};
    else
        use_Flatau = 2; % 1 - using Flatau et al. (1992) coefficients for satruation pressure. 0, 2 - using Alduchov and Eskridge (1996) or Murphy and Koop (2005) empirical equations, respectively, 3 - Teten's formula (ERA-Interim).
    end
else
    use_Flatau = 2; % 1 - using Flatau et al. (1992) coefficients for satruation pressure. 0, 2 - using Alduchov and Eskridge (1996) or Murphy and Koop (2005) empirical equations, respectively, 3 - Teten's formula (ERA-Interim).
end

%% variable determination - FOR TESTS.
% p = 950;
% T = 25;
% RH = 0;

%% Get T in Kelvin
Tk = T + 273.15; %radiosonde temperatures are usually given in Celsius. Therefore, also converting to Kelvin.

%% Calculate latent heat release.
L_v = (2500.8 - 2.36.* T + 0.0016.* T.^2 - 0.00006.* T.^3).* 1000; % latent heat of vaporization in J/Kg (equation is a good approximation for -25 - 40 C range). see Rogers and Yau (1987), Table 2.1.
L_s = (2834.1 - 0.29.* T - 0.004.* T.^2).* 1000; % latent heat of sublimation in J/Kg (almost constant between -40 - 0 C). see Rogers and Yau (1987), Table 2.1.

%% Calculate water vapor saturation pressure
if use_Flatau == 1 % Flatau et al. (1992)
    water_e_s = polyval(p_water_s_coeff, T).*100; % valid for -100 to 0 C (~0.1% max error, T is in Celsius). Based on Flatau et al, J. Applied Met., 1992 (table 5, relative error norm)
    water_e_i = polyval(p_water_i_coeff, T).*100; % valid for -75 to 0 C (~0.384% max error, T is in Celsius). Based on Flatau et al, J. Applied Met., 1992 (table 5, relative error norm)
elseif  use_Flatau == 0 % Alduchov and Eskridge (1996)
    water_e_s = 6.1094 .* exp((17.625 .* T)./ (T + 243.04)).*100; % valid for -40 - +50 C (0.384% max error, T is in Celsius). This equation is empirical based on Alduchov and Eskridge, J. Applied Met., 1996
    water_e_i = 6.1121 .* exp((22.587 .* T)./ (T + 273.16)).*100; % valid for -80 - 0 C (0.414% max error, T is in Celsius). This equation is empirical based on Alduchov and Eskridge, J. Applied Met., 1996
elseif use_Flatau == 2    % Murphy and Koop (2005)
    water_e_s = exp(54.842763 - 6763.22./Tk - 4.210.*log(Tk) + 0.000367.*Tk + tanh(0.0415.*(Tk-218.8)).*(53.878 - 1331.22./Tk - 9.44523.*log(Tk) + 0.014025.*Tk)); % valid for 123 < T < 332 K
    water_e_i = exp(9.550426 - 5723.265./Tk + 3.53068.*log(Tk) - 0.00728332.*Tk); % valid for T > 110 K.
elseif  use_Flatau == 3 % Teten's formula (ERA-Interim)
    water_e_s = 6.1121 .* exp((17.502 .* T)./ (T + 240.97)).*100; % valid for -40 - +50 C (0.384% max error, T is in Celsius). This equation is empirical based on Alduchov and Eskridge, J. Applied Met., 1996
    water_e_i = 6.1121 .* exp((22.587 .* T)./ (T + 273.86)).*100; % valid for -80 - 0 C (0.414% max error, T is in Celsius). This equation is empirical based on Alduchov and Eskridge, J. Applied Met., 1996
elseif  use_Flatau == 4 % Emanuel (1994) (AMPS).
    water_e_s = 6.112 .* exp((17.67 .* T)./ (T + 243.5)).*100; % valid for -40 - +50 C (0.384% max error, T is in Celsius). This equation is empirical based on Alduchov and Eskridge, J. Applied Met., 1996
    water_e_i = exp(23.33086 - 6111.72784./Tk + 0.15215.*log(Tk)).*100; % valid for -80 - 0 C (0.414% max error, T is in Celsius). This equation is empirical based on Alduchov and Eskridge, J. Applied Met., 1996
end

%% Calculate water vapor pressure and mixing ratio.
if exist('water_w', 'var')
    water_e = water_w.*p./ (0.62197 + water_w); % getting p in hPa
    RH = water_e./ water_e_s.* 100;
else
    water_e = RH .* water_e_s ./ 100; % get water vapor pressure from the RH and saturation water pressure. see: http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf.
    water_w = 0.62197 .* water_e./ (p - water_e); % mixing ratio in g/kg. see: http://www.srh.noaa.gov/images/epz/wxcalc/mixingRatio.pdf.
end

%% Calculate potential temperature.
p_d = p - water_e; % Dry air partial pressure.
Tv = Tk.*(1 + (water_w.*1e-3)./0.62197)./(1 + (water_w.*1e-3));
Theta = Tk.* (p_0./ p) .^(R_d/ c_p); % potential temperature.
Theta_e = Tk.* (p_0 ./ p_d).^ (R_d/ c_p) .* (RH./ 100).^ ((water_w).* R_v ./ c_p) .* exp(L_v.* (water_w)./ (c_p.*Tk)); % equivalent potential temperature (neglecting r*c, i.e., total water mixing ratio and heat capacity of liquid water, while retaining good accuracy). see: http://glossary.ametsoc.org/wiki/Equivalent_potential_temperature
Theta_v = Tv.* (p_0./ p) .^(R_d/ c_p); % virtual potential temperature.
%% Calculate first-order equivalent potential temperature eroror.
if exist('delta_Tk', 'var')
    water_e_un = delta_RH .* water_e_s ./ 100;      delta_water_w = 0.62197.* p ./ (water_e - p).^2 .* water_e_un;
    Theta_e_un =  sqrt((((p_0 ./ p_d).^ (R_d/ c_p) .* (RH./ 100).^ ((water_w).* R_v ./ c_p) .* (Tk - L_v.* (water_w)./ c_p).* exp(L_v.* (water_w)./ (c_p.*Tk)) ./ Tk) .* delta_Tk).^2 + ...
        (Tk.* (p_0 ./ p_d).^ (R_d/ c_p) .* (RH./ 100).^ ((water_w).* R_v ./ c_p) .* exp(L_v.* (water_w)./ (c_p.*Tk)).* ((water_w).* R_v ./ c_p)./ RH .* delta_RH).^2 + ...
        (Tk.* (p_0 ./ p_d).^ (R_d/ c_p) .* (RH./ 100).^ ((water_w).* R_v ./ c_p) .* exp(L_v.* (water_w)./ (c_p.*Tk)).* ((L_v./ (c_p.*Tk)) + log(RH./ 100).* (R_v ./ c_p)) .* delta_water_w).^2 );
end
%% Calculate RH_above ice.
RH_i = water_e./ water_e_i .* 100;
RH_l = water_e./ water_e_s .* 100;

%% Arranging ouput in a single struct.
Output_struct.L_v = L_v;
Output_struct.L_s = L_s;
Output_struct.e_s = water_e_s; % Pa
Output_struct.e_i = water_e_i; % Pa
Output_struct.e = water_e; % Pa
Output_struct.w = water_w.*1e3; % g/kg
Output_struct.T_v = Tv; % in K
Output_struct.Theta = Theta;% K
Output_struct.Theta_e = Theta_e;
Output_struct.Theta_v = Theta_v; % in K
if exist('delta_Tk', 'var');    Output_struct.water_e_un = delta_water_w.*1e3;    Output_struct.Theta_e_un = Theta_e_un;      end
Output_struct.RH_i = RH_i; % %
Output_struct.RH_l = RH_l; % %
Output_struct.q = Output_struct.w./1e3./(1+Output_struct.w./1e3).*1e3; % g/kg
Output_struct.q_s = 0.622.*Output_struct.e_s./100./(p./100 - Output_struct.e_s./100).*1e3; % g/kg