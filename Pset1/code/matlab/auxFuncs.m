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
            wages   = data.wages(Iwp);
            wages   = wages/geomean(wages);
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
            tol = 1e-10;
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
                Delta = max(abs(workers-pWorkers));
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
            omega = omega.^(1/par.theta);
            omega(omega>0) = omega(omega>0)/geomean(omega(omega>0));
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
        % Function that solves the fix point for the counterfactual
        % Fixed amenities
        function [out,vnew] = countAgg(par,eq,d_mat,agg,pop)
            % Initial population
            init_pop = sum(eq.num_employment);
            % Unpack endogenous objects
            Q           = eq.home_price;
            wage        = eq.omega;
            floorDist   = eq.floorDist;
            % Fixed objects
            b           = eq.b*agg+eq.B*(~agg);
            land        = eq.land;
            K           = eq.varphi.*land.^(1-par.mu);
            % Unpack objects that will change over the loop
            A = eq.A;
            B = eq.B;
            % Define toelrance and parameter of adjustment
            tol = 1e-2;
            adj = 0.34;
            % Create matrix of commuting costs adjusted by elasticity
            dij = exp(-par.theta*par.kappa*d_mat);
            % Begin fix point routine
            Delta = 10;
            counter = 0;
            flag    = true;
            tic;
            error = [];
            while(Delta>tol)
                % Counter
                counter = counter +1;
                % Big Matrix of latent utilities (adjusted)
                Phi_ij = dij.*((B.*(wage'))./(Q.^(1-par.beta))).^par.theta;
                % Get sum of all utilities
                Phi = sum(Phi_ij,"all");
                % Probabilities (unconditional)
                pi_ij = Phi_ij/Phi;
                % Probabilities (conditional)
                pi_ij_i = dij.*(wage').^par.theta;
                pi_ij_i = pi_ij_i./sum(pi_ij_i,2);
                % Tricky part: Get the total population. 
                % Depends on what you are fixing
                tPop = init_pop*(pop) + Phi*(~pop);
                % Compute workers and residents
                workers  = tPop*sum(pi_ij,1)';
                resident = tPop*sum(pi_ij,2);
                % Production
                Y = A.*(workers.^par.alpha).*(floorDist.*K).^(1-par.alpha);
                % Workers income
                v = sum(pi_ij_i.*(wage'),2);
                % Update endogenous objects
                % Wages
                wage_hat = par.alpha*Y./workers;
                % Price of land
                q    = ((1-par.alpha)*Y+(1-par.beta)*v.*resident)./K;
                idx = (A>0 & B==0);
                q(idx) = (1-par.alpha).*Y(idx)./(floorDist(idx).*K(idx));
                idx = A==0 & B>0;
                q(idx) = ((1-par.beta).*v(idx).*resident(idx))./((1-floorDist(idx)).*K(idx));
                % Floor distribution
                f_hat = ((1-par.alpha)*Y)./(q.*K);
                f_hat(A>0 & B==0) = 1;
                f_hat(A==0 & B>0) = 0;
                % Update B if agglomerations are endogenous
                if(agg)
                    % Get new Omega
                    Omega = exp(-par.rho*d_mat)*(resident./land);
                    % Get new amenities
                    nB =  b.*(Omega.^par.eta);
                    B  = [B,nB]*[adj;1-adj];
                end
                % Calculate differences
                Delta = max(mean(abs([(wage-wage_hat),(Q-q),(floorDist-f_hat)])));
                error = [error;Delta];
                % Update guesses
                new = [[wage,wage_hat];[Q,q];[floorDist,f_hat]]*[adj;1-adj];
                % Unpack new wages
                idx  = 1:par.N;
                wage = new(idx);
                % New housing prices
                idx  = idx + par.N; 
                Q    = new(idx);
                % New floor distribution
                idx  = idx + par.N; 
                floorDist = new(idx);
                % Give an exit
                if(counter>5e4)
                    flag = false;
                    break;
                end
            end
            time = toc;
            % Message for when finish the fix point
            if(flag)
                % Display amount of time and iterations
                fprintf('Fix Point solved. Time elapsed %2.2f seconds in %5i iterations .\n',time,counter);
            else
                figure;
                hold on
                plot(error)
                plot(xlim(),[tol,tol],'--k')
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                %-------------------------------------------------
                error('Fix point not solved. Check for tolerance');
            end
            % Save object in table
            out = eq;
            out.B = B;
            out.floorDist = floorDist;
            out.num_employment = workers;
            out.num_residents  = resident;
            out.omega = wage;
            out.home_price = Q;
            if(agg)
                out.Omega = Omega;
            end
            % Return either population and utility
            if(pop)
                % Save population
                vnew(1) = tPop;
                % Save utility
                vnew(2) = (Phi^(1/par.theta))*par.gammaAdj;    
            else
                % Save population
                vnew(1) = tPop;
                % Save utility
                vnew(2) = 1; 
            end
        end
    
    end
end