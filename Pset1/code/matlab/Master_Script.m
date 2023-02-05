%--------------------------------------------------%
%------- Spatial Economics: Problem Set 1 ---------%
%--------------------------------------------------%

%% Housekeeping and parameters
clear all; clc; close all; 

path.figs = '../../Figures/Single Agent/';
path.tabs = '../../Tables/Single Agent/';

% Only works in Jose's computer. To format figures how I like it. 
try
    global c
    myColors();
end

% Estimate theta
est_theta = false;

% Initialize parameters taken from literature 
par.beta    = 0.76;
par.alpha   = 0.8;
par.mu      = 0.75;
par.theta   = 6.5;
par.kappa   = 0.1/par.theta;
par.rho     = 0.76;
par.eta     = 0.15;

% Add to parameters adjustment from gamma function 
par.gammaAdj = gamma((par.theta-1)/par.theta);

%% Intialize and prepare data

% Read data for locations
data = readtable('../../Data/location_characteristics.csv');
% Read data for commuting costs
flows = readtable('../../Data/commuting_cost.csv');
% Define number of locations 
par.N = size(data,1); 
% Prepare commuting cost 
d_mat = reshape(flows.distance,par.N,par.N)';
% Read shapefile for mapAbsol 
chicago  = readgeotable('../../Data/zillow_neighborhoods/chicago_zillow_neighborhoods.shp');
% Keep only location in Chicago (for which we have data)
chicago = chicago(strcmp(chicago.City,"Chicago"),:);
chicago.RegionID = str2double(chicago.RegionID);
idx = ismember(chicago.RegionID,data.regionid);
%------------------------------------------
fprintf('%4i locations drop because we do not observe prices.\n',sum(~idx));
%------------------------------------------
chicago = chicago(idx,:);
% Add land area to data
land_tab = table(chicago.RegionID,chicago.Shape_Area,'VariableNames',{'regionid','land'});
data     = join(land_tab,data);
% Remove unnecessary stuff
clear land_tab idx

% IDs of locations of interest:
% 1. Hyde Park
% 2. East Hyde Park
% 3. The Loop
id = [269586,403352,269593];
[tf,loc] = ismember(chicago.RegionID,id);
[~,p] = sort(loc(tf));
id = find(tf);
id = id(p);
% Remove unused stuff
clear p tf loc

% Add basline data to Chicago shapefile
chicago = join(chicago,data,"LeftKeys","RegionID","RightKeys","regionid");

%% Invertion of the model and analysis 

% Estimate theta using gravity equation
if(est_theta)
    par = gravityTheta(par,flows);
end

% Solve for preference scale for work location
base_eq = auxFuncs.workScale(par,data,d_mat);

% Invert model 
base_eq = auxFuncs.invModel(par,data,d_mat,base_eq);

%% Create maps

% Create data base to analyze equilibrium objects
chicago_eq = join(chicago,base_eq,"LeftKeys","RegionID","RightKeys","regionid");

% Create a variable for land developed
chicago_eq.land_dev = chicago_eq.varphi.*(chicago_eq.land_base_eq.^(1-par.mu));

% Columns for change
chng = [23,24,28,29,31,32];
names = chicago_eq.Properties.VariableNames;


% Add to table for specific variables
for j= 1:length(chng)
    % Standarize variables for presentation
    chicago_eq.(names{chng(j)})  = normalize(chicago_eq.(names{chng(j)}));
end



% Comercial land distribution
figure; 
gx =geoaxes;
geoplot(chicago_eq,ColorVariable="floorDist");
geobasemap('none')
grid('off')
gx.FontSize = 12;
gx.FontName = 'CMU Serif';
gx.ZoomLevel = 10.84;
gx.LatitudeLabel.FontSize = 15;
gx.LongitudeLabel.FontSize = 15;
cb = colorbar;
cb.FontSize = 12;
clim([-1.5,1.5])
export_fig(strcat(path.figs,'Baseline/theta_dist'),'-pdf','-transparent'); 

% Local Productivity
figure; 
gx =geoaxes;
geoplot(chicago_eq,ColorVariable="A")
geobasemap('none')
grid('off')
gx.FontSize = 12;
gx.FontName = 'CMU Serif';
gx.ZoomLevel = 10.84;
gx.LatitudeLabel.FontSize = 15;
gx.LongitudeLabel.FontSize = 15;
cb = colorbar;
cb.FontSize = 12;
clim([-1.5,1.5])
export_fig(strcat(path.figs,'Baseline/prod'),'-pdf','-transparent'); 


% Local Productivity (Unadjusted)
figure; 
gx = geoaxes;
geoplot(chicago_eq,ColorVariable="Araw")
geobasemap('none')
grid('off')
gx.FontSize = 12;
gx.FontName = 'CMU Serif';
gx.ZoomLevel = 10.84;
gx.LatitudeLabel.FontSize = 15;
gx.LongitudeLabel.FontSize = 15;
cb = colorbar;
cb.FontSize = 12;
clim([-1.5,1.5])
export_fig(strcat(path.figs,'Baseline/prod_raw'),'-pdf','-transparent'); 

% Local Amenities
figure; 
gx = geoaxes;
geoplot(chicago_eq,ColorVariable="B")
geobasemap('none')
grid('off')
gx.FontSize = 12;
gx.FontName = 'CMU Serif';
gx.ZoomLevel = 10.84;
gx.LatitudeLabel.FontSize = 15;
gx.LongitudeLabel.FontSize = 15;
cb = colorbar;
cb.FontSize = 12;
clim([-1.5,1.5])
export_fig(strcat(path.figs,'Baseline/amenities'),'-pdf','-transparent'); 


% Local (Endogenous) Amenities
figure; 
gx = geoaxes;
geoplot(chicago_eq,ColorVariable="b")
geobasemap('none')
grid('off')
gx.FontSize = 12;
gx.FontName = 'CMU Serif';
gx.ZoomLevel = 10.84;
gx.LatitudeLabel.FontSize = 15;
gx.LongitudeLabel.FontSize = 15;
cb = colorbar;
cb.FontSize = 12;
clim([-1.5,1.5])
export_fig(strcat(path.figs,'Baseline/eAmenities'),'-pdf','-transparent'); 


% Density of development
figure; 
gx = geoaxes;
geoplot(chicago_eq,ColorVariable="land_dev")
geobasemap('none')
grid('off')
gx.FontSize = 12;
gx.FontName = 'CMU Serif';
gx.ZoomLevel = 10.84;
gx.LatitudeLabel.FontSize = 15;
gx.LongitudeLabel.FontSize = 15;
cb = colorbar;
cb.FontSize = 12;
clim([-1.5,1.5])
export_fig(strcat(path.figs,'Baseline/land_dev'),'-pdf','-transparent'); 


%% Counterfactual exercise (Fixing population)

% Create full equilibrium object
eq = join(base_eq,data);

% Perturbate parameters
eq.A(id(1)) = eq.A(id(1)) + 0.5*std(eq.A);
eq.b(id(1)) = eq.b(id(1)) + 2*std(eq.b);
eq.B(id(1)) = eq.B(id(1)) + 2*std(eq.B);

% Solve new equilibrium
[eq_new, = auxFuncs.countAgg(par,eq,d_mat,true,true);

% Columns for change
chng = [2,5,10,14,18,22];
names = eq_new.Properties.VariableNames;

names(chng)

% Calculate percentual change
pChng = table2array(eq_new(:,chng))./table2array(eq(:,chng))-1;
% Add to table for specific variables
for j= 1:length(chng)
    eq_new.(names{chng(j)})  = pChng(:,j);
end

% Paste shapefile
eq_new = join(chicago,eq_new,"LeftKeys","RegionID","RightKeys","regionid");


% Wages
figure; 
geoplot(eq_new,ColorVariable="omega")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix Population/wages'),'-pdf','-transparent'); 


% Amenities
figure; 
geoplot(eq_new,ColorVariable="B")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix Population/ammenities'),'-pdf','-transparent'); 

% Floor distribution
figure; 
geoplot(eq_new,ColorVariable="floorDist")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix Population/comFloor'),'-pdf','-transparent'); 

% Employment
figure; 
geoplot(eq_new,ColorVariable="num_employment_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix Population/workers'),'-pdf','-transparent'); 

% Residents
figure; 
geoplot(eq_new,ColorVariable="num_residents_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix Population/residents'),'-pdf','-transparent'); 

% Housing
figure; 
geoplot(eq_new,ColorVariable="home_price_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix Population/housing'),'-pdf','-transparent'); 


%% Counterfactual exercise (Fixing utility)

% Create full equilibrium object
eq = join(base_eq,data);

% Perturbate parameters
eq.A(id(1)) = eq.A(id(1)) + 0.5*std(eq.A);
eq.b(id(1)) = eq.b(id(1)) + 2*std(eq.b);
eq.B(id(1)) = eq.B(id(1)) + 2*std(eq.B);

% Solve new equilibrium
eq_new = auxFuncs.countAgg(par,eq,d_mat,true,false);

% Columns for change
chng = [2,5,10,14,18,22];
names = eq_new.Properties.VariableNames;


% Calculate percentual change
pChng = table2array(eq_new(:,chng))./table2array(eq(:,chng))-1;
% Add to table for specific variables
for j= 1:length(chng)
    eq_new.(names{chng(j)})  = pChng(:,j);
end

% Paste shapefile
eq_new = join(chicago,eq_new,"LeftKeys","RegionID","RightKeys","regionid");


% Wages
figure; 
geoplot(eq_new,ColorVariable="omega")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix utility/wages'),'-pdf','-transparent'); 


% Amenities
figure; 
geoplot(eq_new,ColorVariable="B")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix utility/ammenities'),'-pdf','-transparent'); 

% Floor distribution
figure; 
geoplot(eq_new,ColorVariable="floorDist")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix utility/comFloor'),'-pdf','-transparent'); 

% Employment
figure; 
geoplot(eq_new,ColorVariable="num_employment_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix utility/workers'),'-pdf','-transparent'); 

% Residents
figure; 
geoplot(eq_new,ColorVariable="num_residents_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix utility/residents'),'-pdf','-transparent'); 

% Housing
figure; 
geoplot(eq_new,ColorVariable="home_price_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/Fix utility/housing'),'-pdf','-transparent'); 

