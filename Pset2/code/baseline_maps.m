%%% Plot results from baseline %%%
%%% Baseline is RCP 6.0 %%%


cd("/Users/junwong/Downloads/EGGW_Replication_Package");

data = load("Data/derived/results_forward_Warm_RCP6.0_med.mat");
data_nowarm = load("Data/derived/results_forward_noWarm_RCP6.0.mat");

global H H0 amen_util0 C_vect D_vect Africa_vect ...
     prod0 pop0 pop0_dens pop1_dens pop5_dens u0 ...
     trmult_reduced n beta price_energy0 alpha theta Omega HDI0 ...
     emi0_cell emi0 emi0_ff emi0_no_ff emi_no_ff emi_ff_data_tend chi cost_emi_param ...
     a0 a1 a2 a3 b1 b2 b3 S0 S1 S2 S3 S_preind temp_preind ...
     forc_sens forc_noCO2 c1 c2 d1 d2 temp1 temp2 ...
     temp0_local_long temp0_local Delta_temp1 temp_past Delta_temp ...
     temp0_global scaler_temp theta_amen_scen theta_prod_scen ... 
     lbar lambda gamma1 gamma2 eta mu nu ksi conv_usd ...
     pop_low95 pop_low80 pop_med pop_high80 pop_high95 ...
     pop_low95_hist pop_low80_hist pop_med_hist pop_high80_hist pop_high95_hist ...
     emi_ff_data forc_data length_C_vect length_D_vect length_Africa_vect table_prop ...
     emiCO2_RCP stockCO2_RCP forc_RCP temp_global_RCP ...
     price_clean0_world price_clean0 price_fossil0 zeta_clean0 zeta_fossil0 fossil_share epsilon std_epsilon clean_energy_data tCO2_toe ...
     clean0_cell fossil0_cell map aux_lon aux_lat aux_kron share_agri agri_index ...
     color_darkseagreen color_olive color_darkgreen color_darkcyan color_yellowgreen color_greenyellow ...
     natal_param name_dam_vect name_dam_long_vect cost_CO2_data a_norm ...
     natal_data country_data year_data pop_data max_cumCO2 natal0 natal20 ...
     price_fossil0_world price_fossil0_adj

addpath('Functions')

% Initialize model
initialize(8.5,19500,1.6) % benchmark estimation

% Parameters
T = 400;
HT = repmat(H,1,T);
sigma = 1;

% US file 
out = usaExtract();

%%%%%%%%%%%%%%%%%%%%
%%% Process Data %%%
%%%%%%%%%%%%%%%%%%%%
data20 = load('Data/derived/results20_med.mat');
l0 = data20.l0_model;
realgdp0 = data20.realgdp0_model;

% Extract data
realgdp = data.realgdp_Warm;
amenities = data.amen_Warm.^sigma;
realgdp_nowarm = data_nowarm.realgdp_noWarm;
amenities_nowarm = data_nowarm.amen_noWarm.^sigma;
l = data.l_Warm;
l_nowarm = data_nowarm.l_noWarm;
util = data.u_Warm.^sigma ;
util_nowarm = data_nowarm.u_noWarm.^sigma;

realgdp_matrix = cat(3,realgdp,realgdp_nowarm);
amen_matrix = cat(3,amenities,amenities_nowarm);
pop_matrix = cat(3,l,l_nowarm);
util_matrix = cat(3,util,util_nowarm);

% Initialize matrices
realgdp_grid = NaN(180,360,T+1,2);
realgdp_aux = realgdp_grid(:,:,1);
realgdp_aux(H0>0) = realgdp0;
realgdp_grid(:,:,1,:) = repmat(realgdp_aux,1,1,1,2);
amen_grid = NaN(180,360,T+1,2);
amen_aux = amen_grid(:,:,1);
amen_aux(H0>0) = a_norm;
amen_grid(:,:,1,:) = repmat(amen_aux,1,1,1,2);
pop_grid = NaN(180,360,T+1,2);
pop_aux = pop_grid(:,:,1);
pop_aux(H0>0) = pop0_sh;
pop_grid(:,:,1,:) = repmat(pop_aux,1,1,1,2);
util_grid = NaN(180,360,T+1,2);
util_aux = util_grid(:,:,1);
util_aux(H0>0) = u0; 
util_grid(:,:,1,:) = repmat(util_aux,1,1,1,2);

pop_sh = cat(3,l.*HT./repmat(sum(l.*HT),n,1),l_nowarm.*HT./repmat(sum(l_nowarm.*HT),n,1));

% Arrange data into 180 x 360 grid by each year 
for t=2:T+1 
    for i=1:2
        
    realgdp_aux = realgdp_grid(:,:,t,i);
    realgdp_aux(H0>0) = realgdp_matrix(:,t-1,i);
    realgdp_grid(:,:,t,i) = realgdp_aux;
 
    %realgdp_ratio_aux = realgdp_matrix(:,t-1,1)./realgdp_matrix(:,t-1,2);
    %[hist_realgdp(:,t), interv_realgdp(:,t)] = ...
    %    histwc(realgdp_matrix(:,t-1,1)./realgdp_matrix(:,t-1,2), pop_sh(:,t-1), nbins);
	%std_realgdp(t) = sqrt(sum(pop_sh(:,t-1).*((realgdp_ratio_aux-sum(realgdp_ratio_aux.*pop_sh(:,t-1))).^2)));

    amen_aux = amen_grid(:,:,t,i);
    amen_aux(H0>0) = amen_matrix(:,t-1,i);
    amen_grid(:,:,t,i) = amen_aux;

    %amen_ratio_aux = amen_matrix(:,t-1,1)./amen_matrix(:,t-1,2);
    %std_amen(t) = sqrt(sum(pop_sh(:,t-1).*((amen_ratio_aux-sum(amen_ratio_aux.*pop_sh(:,t-1))).^2)));

    pop_aux = pop_grid(:,:,t,i);
    pop_aux(H0>0) = pop_matrix(:,t-1,i);
    pop_grid(:,:,t,i) = pop_aux;

    util_aux = util_grid(:,:,t,i);
    util_aux(H0>0) = util_matrix(:,t-1,i);
    util_grid(:,:,t,i) = util_aux;
    end
end

% Get PDV of real GDP and welfare 
PDV_realgdp_matrix = NaN(n,2);
PDV_realgdp_grid = NaN(180,360,2);
PDV_util_matrix = NaN(n,2);
PDV_util_grid = NaN(180,360,2);
beta_matrix = beta.^repmat(0:1:T-1,n,1);

for i=1:2
    PDV_realgdp_matrix(:,i) = sum(beta_matrix.*realgdp_matrix(:,:,i),2);
    realgdp_w = sum(pop_sh(:,:,i).*realgdp_matrix(:,:,i));
    BGP_realgdp_w(i) = realgdp_w(T)./realgdp_w(T-1);
    index_last_PDV_realgdp = (beta*BGP_realgdp_w(i) < 1);    
    arg_last_PDV_realgdp = (beta^T*realgdp_matrix(:,T,i).*BGP_realgdp_w(i)) ./...
            (1-beta*BGP_realgdp_w(i));
    PDV_realgdp_matrix(:,i) = PDV_realgdp_matrix(:,i) + arg_last_PDV_realgdp;
    if index_last_PDV_realgdp == 0
        PDV_realgdp_matrix(:,i) = NaN;
    end
    PDV_realgdp_grid_aux = NaN(180,360); 
    PDV_realgdp_grid_aux(H0>0) = PDV_realgdp_matrix(:,i);
    PDV_realgdp_grid(:,:,i) = PDV_realgdp_grid_aux;

    PDV_util_matrix(:,i) = sum(beta_matrix.*util_matrix(:,:,i),2);
    util_w = sum(pop_sh(:,:,i).*util_matrix(:,:,i));
    BGP_util_w(i) = util_w(T)./util_w(T-1);
    index_last_PDV_util = (beta*BGP_util_w(i) < 1);    
    arg_last_PDV_util = (beta^T*util_matrix(:,T,i).*BGP_util_w(i)) ./...
            (1-beta*BGP_util_w(i));
    PDV_util_matrix(:,i) = PDV_util_matrix(:,i) + arg_last_PDV_util;
    if index_last_PDV_util == 0
        PDV_util_matrix(:,i) = NaN;
    end
    PDV_util_grid_aux = NaN(180,360); 
    PDV_util_grid_aux(H0>0) = PDV_util_matrix(:,i);
    PDV_util_grid(:,:,i) = PDV_util_grid_aux;
end

%%%%%%%%%%%%%%%%%%%%%%
%%% Plot amenities %%%
%%%%%%%%%%%%%%%%%%%%%%
for t=[51 101 201]

    amen_data = map.*kron(amen_grid(:,:,t,1)./amen_grid(:,:,t,2),aux_kron);
    us_coords = map.*kron(out, aux_kron);
    
    amen_data(us_coords==0) = NaN;
    [r, v] = find(~isnan(amen_data));
    amen_data = amen_data(min(r):max(r), min(v):1800);

    figure
    colormap(jet);
	imagesc(aux_lon,aux_lat,...
    	amen_data,...
        'AlphaData',~isnan(amen_data), [0.94 1.1]); 
    colorbar('eastoutside');
    set(gcf,'color','w');
    set(gca,'color',[1 1 1],'FontSize',14)
    %set(gca,'Xtick',[-150 -100 -50 0 50 100 150],'XTickLabel',[-150 -100 -50 0 50 100 150])    
    %set(gca,'Ytick',[-80 -60 -40 -20 0 20 40 60 80],'YTickLabel',[-80 -60 -40 -20 0 20 40 60 80])    
    title(['Amenities in ' num2str(1999+t) ': RCP 6.0' ' relative to no warming'])
    ax = gca;
    ax.YTickLabel = flipud(ax.YTickLabel);
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3) - 0.025;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    saveas(gcf,['Maps and Figures/amen_' num2str(1999+t) '.pdf'])
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot population %%%
%%%%%%%%%%%%%%%%%%%%%%%
pop0_sh = pop0_dens.*H/lbar;

for t=[51 101 201]

    pop_data = map.*kron(pop_grid(:,:,t,1)./pop_grid(:,:,t,2),aux_kron);
    us_coords = map.*kron(out, aux_kron);
    
    pop_data(us_coords==0) = NaN;
    [r, v] = find(~isnan(pop_data));
    pop_data = pop_data(min(r):max(r), min(v):1800);

    figure
    colormap(jet);
	imagesc(aux_lon,aux_lat,...
    	pop_data,...
        'AlphaData',~isnan(pop_data), [0.94, 1.2]); 
    colorbar('eastoutside');
    set(gcf,'color','w');
    set(gca,'color',[1 1 1],'FontSize',14)
    %set(gca,'Xtick',[-150 -100 -50 0 50 100 150],'XTickLabel',[-150 -100 -50 0 50 100 150])    
    %set(gca,'Ytick',[-80 -60 -40 -20 0 20 40 60 80],'YTickLabel',[-80 -60 -40 -20 0 20 40 60 80])    
    title(['Population Density in ' num2str(1999+t) ': RCP 6.0' ' relative to no warming'])
    ax = gca;
    ax.YTickLabel = flipud(ax.YTickLabel);
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3) - 0.025;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    saveas(gcf,['Maps and Figures/pop_' num2str(1999+t) '_' '.pdf'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot PDV real GDP %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

realgdp_data = map.*kron(PDV_realgdp_grid(:,:,1)./PDV_realgdp_grid(:,:,2),aux_kron);
us_coords = map.*kron(out, aux_kron);
    
realgdp_data(us_coords==0) = NaN;
[r, v] = find(~isnan(realgdp_data));
realgdp_data = realgdp_data(min(r):max(r), min(v):1800);

figure
colormap(jet);
imagesc(aux_lon,aux_lat,realgdp_data,...
	'AlphaData',~isnan(realgdp_data)); 
colorbar('eastoutside');
saveas(gcf,'Maps and Figures/pdv_gdp.pdf')

%%%%%%%%%%%%%%%%%%%%
%%% Plot welfare %%%
%%%%%%%%%%%%%%%%%%%%
util_data = map.*kron(PDV_util_grid(:,:,1)./PDV_util_grid(:,:,2),aux_kron);
us_coords = map.*kron(out, aux_kron);
    
util_data(us_coords==0) = NaN;
[r, v] = find(~isnan(util_data));
util_data = util_data(min(r):max(r), min(v):1800);

figure
colormap(jet);
imagesc(aux_lon,aux_lat,util_data,...
	'AlphaData',~isnan(util_data)); 
colorbar('eastoutside');

saveas(gcf,'Maps and Figures/pdv_welfare.pdf')
