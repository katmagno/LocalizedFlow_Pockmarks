close all;
clear all;
% paths = [{'/Users/kmagno/Documents/porosity_update/np21/'},{'/Users/kmagno/Documents/porosity_update/np24/'},{'/Users/kmagno/Documents/porosity_update/np27/'},{'/Users/kmagno/Documents/porosity_update/np30/'},...
%     {'/Users/kmagno/Documents/porosity_update/np33/'}];
paths = [DATA PATHS];

% Physical Parameters
perm = 5e-15; % m^2
mu_f = 5e-5; % viscosity of fluid, Pas
eta_phi = 5e15; % Bulk viscosity, Pas
delta_rhog = 1000*9.81; % Buoyancy phase differential
tau_c = ((eta_phi*mu_f/perm)^0.5)/delta_rhog; % Characteristic time
delta_c = sqrt(perm*eta_phi/mu_f); % Characteristic length
tau_years = tau_c/(3600*24*365); % Characteristic time
flux_scale = (delta_c/tau_years); % Characteristic flux

% Allocate space
vol_flux = zeros([13,5]);
error_volflux = zeros([13,5]);
% Loop through permeability power-law exponents
for jj = 1:length(paths)
    disp('Path =')
    disp(jj)
    clearvars -except mean_test delta_c tau_years flux_scale paths jj kk perm data_count vol_flux error_volflux delta_ny cr_flux cr_area
    % Set model resolution
    res = 256;
    ny = res*2-1;
    nx = res-1;
    ly = 70;
    lx = ly/2;
    dy = ly/(ny) ;
    % Convert 200m to model resolution units
    delta_ny = (200/delta_c)*(ny/ly);

    syms A x y x0 y0 sigmaX sigmaY a b c d
    f(A,x,y,x0,y0,sigmaX,sigmaY, a, b, c, d) = A * exp(-((x-x0).^2./(2*sigmaX^2)+(y-y0).^2./(2*sigmaY^2)));
    I1 = int(f,x, a, b);
    I2 = int(I1,y, c, d); % Indefinite 2D integral
    fA(A,x0,y0,sigmaX,sigmaY) = diff(I2,A);
    fx0(A,x0,y0,sigmaX,sigmaY) = diff(I2,x0);
    fy0(A,x0,y0,sigmaX,sigmaY) = diff(I2,y0);
    fsigmaX(A,x0,y0,sigmaX,sigmaY) = diff(I2,sigmaX);
    fsigmaY(A,x0,y0,sigmaX,sigmaY) = diff(I2,sigmaY);

    % Set plot colors and font size
    colors_256_3 = [191 211 230; 158 188 218; 140 150 198; 136 86 167; 129 15 124];
    colors_256_3_darker = ([191 211 230; 158 188 218; 140 150 198; 136 86 167; 129 15 124]).*0.90;
    colors_3= colors_256_3./256;
    colors_3_darker = colors_256_3_darker./256;
    FS = 18;

    % Import data
    import_path = string(paths{jj});
    color_num = jj;
    time_cell = dir(strcat(import_path,'times*.csv'));
    time_fnm = {time_cell(:).name};
    time_fnm_oi = string(time_fnm{1});
    times_t = readmatrix(strcat(import_path,time_fnm_oi));
    A = dir(strcat(import_path,'qDys*.csv'));
    names = {A(:).name};
    B=natsort(names);
    adj_t = times_t(end)/length(B):times_t(end)/length(B):times_t(end);
    P = dir(strcat(import_path,'Phi*.csv'));
    names_phi = {P(:).name};
    Phi_fnm=natsort(names_phi);

    count = 0;
    for ii = 1:length(B)
        T = readmatrix(strcat(import_path,string(B{ii})));
        T = rot90(T);
        P = readmatrix(strcat(import_path,string(Phi_fnm{ii})));
        Phi = rot90(P);
        % Identify location of seafloor, horz_slice
        if ii == 1
            [rs,~] = find(mean(T,2) > 1.0);
            horz_slice = round(rs(1)-delta_ny);
        end
        slices(ii,:) = T(horz_slice,:);
        phi_slices(ii,:) = Phi(horz_slice,:);
    end
    % Threshold flux data to ignore fluxes lower than 2.
    slices(slices < 2)=0;
    % Find where the pockmarks form
    K = find(sum(slices,2) ~=0);
    zero_ind = find(slices==0);
    start_oi = K(1);
    rows_oi = slices(start_oi:end,:);
    figure(12); imagesc(rows_oi);
    
    times_oi = adj_t(start_oi:end);
    phi_oi = phi_slices(start_oi:end,:);
    time_period = times_oi.*tau_years;
    tot_real_time = time_period(end);

    locs_list = [];
    full_pic = zeros(2,2);
    full_rad = zeros(2,2);
    flux_std = zeros(2,2);
    mean_phi = zeros(2,2);
    area_std = zeros(2,2);
    full_vol = zeros(2,2);
    vol_err = zeros(2,2);
    full_phi = zeros(2,2);
    count_loc = 0;
    % find pockmarks in flux data along seafloor
    time_0 = 0;
    for rr = 1:size(rows_oi,1)
        disp('rr = ')
        disp(round(rr/size(rows_oi,1),3))
        dt = times_oi(rr)-time_0;
        single_row = rows_oi(rr,:);
        single_phi = (phi_oi(rr,:));
        [pks,locs,w,proms] = findpeaks(single_row,'WidthReference','halfheight');
        for nn = 1:length(locs)
            %Determine area of each peak for each loc
            if proms(nn) > 1.5
                start_pt = floor(locs(nn)-(w(nn)/2));
                end_pt = floor(locs(nn)+(w(nn)/2));
                if start_pt <= 0
                    start_pt = 1;
                end
                if end_pt > nx
                    end_pt = nx;
                end
                int_row = single_row(start_pt:end_pt);
                int_phi = single_phi(start_pt:end_pt);
                if length(find(int_row ~=0)) < 3
                    continue
                else
                  
                    int_flux = mean(int_row);
                    mean_phi = mean(single_phi);
                    pk_radius = w(nn)/2;
                    v = (1:length(int_row))*delta_c; % m
                    vnd = (1:length(int_row)); 
                    gauss_func = fit(v.',(int_row*flux_scale).','gauss1');
                    ci = confint(gauss_func,0.95);
                    [vol, dI] = calc_vol(v, gauss_func, ci, fA, fx0, fy0, fsigmaX, fsigmaY);
                    gauss_vol = vol*dt*tau_years;
                    err_vol = dI*dt*tau_years;
                    FWHM = 2 * sqrt(2*log(2)) * gauss_func.c1;
                    if pk_radius == 0
                        continue
                    else
                        if ~any(locs(nn) == locs_list)
                            if (any(locs(nn) == locs(1:nn-1)+1) || any(locs(nn) == locs(1:nn-1)-1))
                                continue
                            else
                                count_loc = count_loc + 1;
                                locs_list(count_loc) = locs(nn);
                                full_vol(1,count_loc) = gauss_vol;
                                vol_err(1,count_loc) = err_vol;
                                full_pic(1,count_loc) = int_flux;
                                flux_std(1,count_loc) = std(int_row);
                                full_rad(1,count_loc) = pk_radius;
                                area_std(1,count_loc) = 1/2;
                                full_phi(1,count_loc) = mean_phi;
                            end
                        elseif any(locs(nn) == locs_list)
                            count_loc_copy = find(locs(nn)==locs_list);
                            add_row = size(full_pic(:,count_loc_copy),1);
                            full_pic(add_row+1,count_loc_copy) = int_flux;
                            full_vol(add_row+1,count_loc_copy) = gauss_vol;
                            vol_err(add_row+1,count_loc_copy) = err_vol;
                            flux_std(add_row+1,count_loc_copy) = std(int_row);
                            full_rad(add_row+1,count_loc_copy) = pk_radius;
                            area_std(add_row+1,count_loc_copy) = 1/2;
                            full_phi(add_row+1,count_loc_copy)= mean_phi;
                        end
                    end
                end
            else
                continue
            end
        end
        time_0 = times_oi(rr);
    end
    % Convert to radius to area for pockmarks, assume circle
    count_std = 0;
    for gg = 1:size(full_vol,2)
        if nnz(nonzeros(full_vol(:,gg))) == 0
            continue
        else
            nz_avg_flux(gg) = mean(nonzeros(full_pic(:,gg)));
            nz_tot_vol(gg) = sum(nonzeros(full_vol(:,gg))); % totl vol per pockmark
            nz_var_vol(gg) = sum(nonzeros(vol_err(:,gg)));
            nz_avg_radius(gg) = mean(nonzeros(full_rad(:,gg)));
            nz_var_rad(gg) = 0.5./sqrt(length(nonzeros(area_std(:,gg))));
            nz_var_flux(gg) = sqrt(sum((nonzeros(flux_std(:,gg))).^2)/length(nonzeros(flux_std(:,gg))));
            nz_phi(gg) = mean(nonzeros(full_phi(:,gg)));
        end
    end

    SEM_rad= nz_var_rad;
    SEM_flux= nz_var_flux;
    fin_flux= nz_avg_flux;
    fin_rad = nz_avg_radius;
    fin_phi = nz_phi;

    fin_flux_kmpyr = fin_flux*flux_scale*1e-3; %km/yr
    fin_Sigmaflux = SEM_flux*flux_scale*1e-3; %km/yr
    fin_vol = nz_tot_vol*1e-9; % km^3
    fin_vol_err = nz_var_vol*1e-9;


    fin_SigmaArea = (pi*(SEM_rad*(lx/nx)*delta_c).^2)*1e-6;
    fin_area_km2 = (pi*(fin_rad*(lx/nx)*delta_c).^2)*1e-6; %km^2

    figure(6),set(gcf,'Position',[45 99 1066 663]);
    err1_color = [254 224 182]./255;
    err2_color = [253 184 99]./255;
    scatter(fin_area_km2,fin_flux_kmpyr,150*ones([1 length(fin_flux_kmpyr)]),'filled','LineWidth',60,'MarkerFaceColor',colors_3(jj,:),'MarkerEdgeColor',colors_3_darker(jj,:),'LineWidth',2); grid on;
    hold on;
    leg= legend('$\textbf{\boldmath{n$_{k}$=2.1}}$',...
        '$\textbf{\boldmath{n$_{k}$=2.4}}$',...
        '$\textbf{\boldmath{n$_{k}$=2.7}}$',...
        '$\textbf{\boldmath{n$_{k}$=3.0}}$','$\textbf{\boldmath{n$_{k}$=3.3}}$'); ...
    set(leg,'interpreter','latex','Fontsize',FS-3,'FontWeight','bold');
    ax2 = gca;
    ax2.FontSize = FS;
    ax2.TickLabelInterpreter = 'latex';
    xlabel('$\textbf{Pockmark Area \boldmath{[km$^{2}$]}}$','interpreter','latex','Fontsize',FS);
    ylabel('$\textbf{Fluid Release Rate \boldmath{[km/yr]}}$','interpreter','latex','Fontsize',FS,'FontWeight','bold');
    hold on; 
    % Sort sizes of simulated pockmarks according to observed pockmark
    % sizes on the Chatham Rise
     p_flux = 0;
     p_flux_var =0;
     p_area = 0;
     p_area_var = 0;
     p_phi=0;
     p_vol_rate = 0;
     p_vol_err = 0;

    for ff = 1:size(fin_area_km2,2)
        area_oi = fin_area_km2(ff);
        if area_oi < 0.05
            continue
        elseif area_oi > 1.45
            continue
        elseif 0.05 < round(area_oi,2) && round(area_oi,2) <= 0.15
            aa = 1;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.15 < round(area_oi,2) && round(area_oi,2) <= 0.25
            aa = 2;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.25 < round(area_oi,2) && round(area_oi,2)<= 0.35
            aa = 3;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.35 < round(area_oi,2) && round(area_oi,2) <= 0.45
            aa = 4;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.45 < round(area_oi,2) && round(area_oi,2) <= 0.55
            aa = 5;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.55 < round(area_oi,2) && round(area_oi,2) <= 0.65
            aa = 6;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.65 < round(area_oi,2) && round(area_oi,2) <= 0.75
            aa = 7;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.75 < round(area_oi,2) && round(area_oi,2) <= 0.85
            aa = 8;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 0.85 < round(area_oi,2) && round(area_oi,2) <= 0.95
            aa = 9;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 1.05 < round(area_oi,2) && round(area_oi,2) <= 1.15
            aa = 10;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 1.15 < round(area_oi,2) && round(area_oi,2) <=  1.25
            aa = 11;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 1.25 < round(area_oi,2) && round(area_oi,2) <= 1.35
            aa = 12;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        elseif 1.35 < round(area_oi,2) && round(area_oi,2) <= 1.45
            aa = 13;
            p_flux(aa,end+1) = fin_flux_kmpyr(ff);
            p_flux_var(aa,end+1) = fin_Sigmaflux(ff);
            p_phi(aa,end+1) = fin_phi(ff);
            p_area(aa,end+1)  = fin_area_km2(ff);
            p_area_var(aa,end+1) = fin_SigmaArea(ff);
            p_vol_rate(aa,end+1) = fin_vol(ff);
            p_vol_err(aa,end+1) = fin_vol_err(ff);

        end
    end
    % Determine fluid volume transfer rates scaled by porosity
    for gg = 1:size(p_vol_rate,1) %size class x pockmark #
        if nnz(nonzeros(p_vol_rate(gg,:))) == 0
            vol_flux(gg,jj) = 0;
        else
        vol_flux(gg,jj) = mean(nonzeros(p_vol_rate(gg,:)));
        error_volflux(gg,jj) = sqrt(sum(p_vol_err(gg,:).^2)) / sqrt(length(p_vol_err(gg,:)));
        cr_flux(gg,jj) = mean(nonzeros(p_flux(gg,:)));%average over pockmark numbers for each size class
        cr_area(gg,jj) = mean(nonzeros(p_area(gg,:)));
        end
    end

end

%% Scale fluid volume transfer rates by pockmark size distributions on CR
obs_pockmark_n = [40; 114; 136; 89; 52; 23; 11; 1; 4; 2; 1; 1; 2];
cr_vol_rate = vol_flux.*obs_pockmark_n;
cr_vol_err = error_volflux .*obs_pockmark_n;
final_vol_flux = mean(sum(cr_vol_rate,1))/tot_real_time %km3/yr
sig_sum = sqrt(sum(cr_vol_err.^2,1));
sig_mean = sqrt(sum(sig_sum.^2))/length(sig_sum);
final_vol_sig = sig_mean/tot_real_time
% Convert to mass in PgCO2
% liquid CO2 kg/km^3
rho_l = 1.1e12; 
% gas CO2
rho_g = 1.87e+9;

kg2Pg = 1e-12; % Convert kg to Pg

mass_flux_liq = final_vol_flux*rho_l*kg2Pg;
mass_flux_gas = final_vol_flux*rho_g*kg2Pg;

disp('Mass Flux for Liquid:')
disp(mass_flux_liq)

disp('Mass Flux for Gas:')
disp(mass_flux_gas)

disp('DONE')