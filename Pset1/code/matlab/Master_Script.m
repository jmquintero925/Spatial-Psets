%--------------------------------------------------%
%------- Spatial Economics: Problem Set 1 ---------%
%--------------------------------------------------%

%% Housekeeping and parameters
clear all; clc; close all; 

path.figs = '../../Figures/';
path.tabs = '../../Tables/';

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


% Comercial land distribution
figure; 
geoplot(chicago_eq,ColorVariable="floorDist")
geobasemap('none')
grid('off')
cb = colorbar;
set(gca,'ColorScale','log')
export_fig(strcat(path.figs,'theta_dist'),'-pdf','-transparent'); 

% Local Productivity
figure; 
geoplot(chicago_eq,ColorVariable="A")
geobasemap('none')
grid('off')
cb = colorbar;
set(gca,'ColorScale','log')
export_fig(strcat(path.figs,'prod'),'-pdf','-transparent'); 


% Local Productivity (Unadjusted)
figure; 
geoplot(chicago_eq,ColorVariable="Araw")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'prod_raw'),'-pdf','-transparent'); 

% Local Amenities
figure; 
geoplot(chicago_eq,ColorVariable="B")
geobasemap('none')
grid('off')
cb = colorbar;
set(gca,'ColorScale','log')
export_fig(strcat(path.figs,'amenities'),'-pdf','-transparent'); 


% Local (Endogenous) Amenities
figure; 
geoplot(chicago_eq,ColorVariable="B")
geobasemap('none')
grid('off')
cb = colorbar;
set(gca,'ColorScale','log')
export_fig(strcat(path.figs,'eAmenities'),'-pdf','-transparent'); 


% Density of development
figure; 
geoplot(chicago_eq,ColorVariable="varphi")
geobasemap('none')
grid('off')
cb = colorbar;
set(gca,'ColorScale','log')
export_fig(strcat(path.figs,'varphi'),'-pdf','-transparent'); 


%% Counterfactual exercise

% Create full equilibrium object
eq = join(base_eq,data);

% Perturbate parameters
eq.A(id(1)) = eq.A(id(1))*2;
eq.b(id(1)) = eq.b(id(1))*2;
eq.B(id(1)) = eq.B(id(1))*2;

% Solve new equilibrium
eq_new = auxFuncs.countAgg(par,eq,d_mat,true,true);

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


%% Create maps

% Wages
figure; 
geoplot(eq_new,ColorVariable="omega")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/wages'),'-pdf','-transparent'); 


% Amenities
figure; 
geoplot(eq_new,ColorVariable="B")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/ammenities'),'-pdf','-transparent'); 

% Floor distribution
figure; 
geoplot(eq_new,ColorVariable="floorDist")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/comFloor'),'-pdf','-transparent'); 

% Employment
figure; 
geoplot(eq_new,ColorVariable="num_employment_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/workers'),'-pdf','-transparent'); 

% Residents
figure; 
geoplot(eq_new,ColorVariable="num_residents_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/residents'),'-pdf','-transparent'); 

% Housing
figure; 
geoplot(eq_new,ColorVariable="home_price_eq_new")
geobasemap('none')
grid('off')
cb = colorbar;
export_fig(strcat(path.figs,'Counterfactual/housing'),'-pdf','-transparent'); 


