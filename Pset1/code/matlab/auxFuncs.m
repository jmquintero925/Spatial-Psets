% File wirh auxiliary functions for the
classdef    auxFuncs
    methods     ( Static = true )
        % Function that inverts the model
        function eq_tab = invModel(par,data,d_mat,eq_tab)
            tic;
            % Locations with workers and total locations with workers
            Iwp = (data.num_employment>0);
            % Locations with residents
            Irp = (data.num_residents>0);
            % Get workers
            workers     = data.num_employment(Iwp);
            % Get residents
            resident    = data.num_residents(Irp);
            % Get wages 
            wages   = eq_tab.omega(Iwp);
            % Get commuting costs
            d_ij = exp(-par.theta*par.kappa*d_mat(Irp,Iwp));
            % Normalize residents relative to geometric mean
            residents_n = resident/sum(resident);
            residents_n = residents_n/geomean(residents_n);
            % Floor prices normalized (In locations with residents)
            q_n = data.home_price(Irp)/geomean(data.home_price(Irp));
            % Commuting access normalized
            com_acces = d_ij*(wages.^par.theta);
            com_acces = com_acces/geomean(com_acces);
            % Amenities
            B   = (residents_n.^(1/par.theta)).*(q_n.^(1-par.beta))./(com_acces.^(1/par.theta));
            %----- Add to structure ------%
            eq_tab.B = zeros(par.N,1);
            eq_tab.B(Irp) = B;
            % Productivity
            A   = ((data.home_price(Iwp)/(1-par.alpha)).^(1-par.alpha)).*((wages/par.alpha).^par.alpha);
            %----- Add to structure ------%
            eq_tab.A      = nan(par.N,1);
            eq_tab.A(Iwp) = A; 
            % Conditional probabilities of living in i
            pi_iji = (d_ij.*(wages').^par.theta);
            pi_iji = pi_iji./sum(pi_iji,2);
            % Land residential
            L_r     = (1-par.beta)*(pi_iji*wages).*resident./data.home_price(Irp);
            %----- Add to structure ------%
            eq_tab.Lr   = zeros(par.N,1);
            eq_tab.Lr(Irp) = L_r;
            % Land production
            L_m = ((wages./(par.alpha*A)).^(1/(1-par.alpha))).*workers;
            %----- Add to structure ------%
            eq_tab.Lm   = zeros(par.N,1);
            eq_tab.Lm(Irp) = L_m;
            % Density of development
            eq_tab.varphi = (eq_tab.Lm + eq_tab.Lr).*(data.land.^(par.mu-1));
            % Distribution of floor
            eq_tab.floorDist = (eq_tab.Lm)./(eq_tab.Lm+eq_tab.Lr);
            % Raw productivity 
            eq_tab.Araw      = eq_tab.A./(eq_tab.E.^(par.alpha/par.theta)); 
            % Spillovers of externalities
            Omega  = exp(-par.rho*d_mat(Irp,Irp))*(resident./data.land(Irp));
            %----- Add to structure ------%
            eq_tab.Omega(Irp) = Omega; 
            % Endogenous Amenities
            eq_tab.b    = B./(Omega.^par.eta);
            time = toc;
            fprintf('Model Solved. Time elapsed %2.2f.\n',time);
        end
        % Function that calculates the scale parameter for work location
        function eq_tab = workScale(par,data,d_mat)
            eq_tab = table();
            % Initialize vector of preferences
            eq_tab.E        = zeros(par.N,1);
            eq_tab.omega    = zeros(par.N,1);
            % Locations with workers and total locations with workers
            Iwp         = (data.num_employment>0);
            workers     = data.num_employment(Iwp);
            residents   = data.num_residents;
            % Initialize adjusted vector of adjusted wages wages
            omega = ones(sum(Iwp),1);
            % Define toelrance and parameter of adjustment
            tol = 1e-06;
            adj = 0.5;
            % Create matrix of commuting costs adjusted by elasticity
            dij = exp(-par.theta*par.kappa*d_mat(:,Iwp));
            % Begin fix point routine
            Delta = 10;
            counter = 0;
            flag    = true;
            tic;
            while(Delta>tol)
                % Add iteration
                counter = counter+1;
                % Construct matrix
                C_omega = (omega').*dij;
                C_omega = C_omega./sum(C_omega,2);
                % Get predicted values of residents
                pWorkers = (C_omega')*residents;
                % Get new value of Delta
                Delta = max((workers-pWorkers).^2);
                % Update wages (adjusted) if necessary
                omega_hat = (workers./pWorkers).*omega;
                omega = [omega,omega_hat]*[adj;1-adj];
                if(counter==5e3)
                    flag = false;
                    break;
                end
            end
            time = toc;
            % Message for when finish the fix point
            if(flag)
                % Display amount of time and iterations
                fprintf('Fix Point solved. Time elapsed %2.2f,\t Number of iterations %5i.\n',time,counter);
            else
                error('Fix point not solved. Check for tolerance');
            end
            % Back up preference parameters
            E = (omega./data.wages(Iwp)).^(par.theta);
            % Save in structure
            eq_tab.E(Iwp)     = E;
            eq_tab.omega(Iwp) = omega;
            % Add other objects to equilibrium object
            eq_tab.regionid = data.regionid;
            eq_tab.land = data.land;
        end
        % Gravity equation to etimate theta from commuting flows
        function par = gravityTheta(par,flows)
            
        end
    
    end
end