
function [Best, FEs] = CLDMSL_PSO(Info, para, CostFunction)
rand('state',sum(100*clock));

fhd1=str2func('cec19_func');

% -------------- Extract Information ----------------------------
func_num            = Info.func_num;
max_FEs             = Info.max_FEs;
dimension           = Info.dimension;
f_star              = Info.f_star;

% -------------- The domain of the function ---------------------
Xmin                = Info.Search_Range(1);
Xmax                = Info.Search_Range(2);

% -------------- Initial Parameters of CLDMSL_PSO ---------------------
num_particle        = para.PopSize;         
max_iteration = ceil(max_FEs/num_particle);     % Max. no. of iteration number

% *********************** Parameter set ***********************************
num_g = num_particle;
num_g1 = 20;
num_g2 = num_g - num_g1;           
group_ps = 5;                  
group_num = num_g2/group_ps;    

j = 0:(1/(num_g - 1)):1;          
j = j*10;
Pc = ones(dimension,1)*(0.05 + ((0.45).*(exp(j) - exp(j(1)))./(exp(j(num_g)) - exp(j(1)))));      
Weight = 0.99 - (1:max_iteration)*0.79/max_iteration;                
Weight1 = 0.2 + (0.99-0.2)*(1./(1 + exp(-5*(2*(max_iteration - (1:max_iteration))/max_iteration - 1)))); 
C = 0.1;                                                                            
c = 3 - (1:max_iteration)*1.5/max_iteration;                           
c1 = 2.5 - (1:max_iteration)*2/max_iteration;                  
c2 = 0.5 + (1:max_iteration)*2/max_iteration; 
R = ceil(10+(30-10).*((max_iteration-(1:max_iteration))./max_iteration));

flag = 0;  gap1 = 5; sigmax = 1;   sigmin = 0.1;  sig = 1;   

L_FEs = 100;
L_num = ceil(0.25.*group_num);
cc=[1.49445 1.49445]; 
iwt=0.7298;
LB = Xmin;
UB = Xmax;

% **********************Initialization*************************************
if length(Xmin) == 1
    Xmin = repmat(Xmin,num_g,dimension);       
    Xmax = repmat(Xmax,num_g,dimension);
else
    Xmin = repmat(Xmin,num_g,1);
    Xmax = repmat(Xmax,num_g,1);
end
v_max = 0.2*(Xmax-Xmin);                         
v_min = -v_max;

%% initialize the position and velocity of the particles
pos = Xmin + (Xmax-Xmin).*rand(num_g,dimension);     % position       
vel = v_min + 2.*v_max.*rand(num_g,dimension);     % velocity 

k = 0;     FEs = 0;
for i = 1:num_g  
       result(i,1) = CostFunction(pos(i,:)',func_num); % calculation fitness function
end
FEs = num_g;

%% initialize the pbest and the pbest's fitness value
pbest_pos = pos;                            
pbest_val = result';      
[gbest_val,g_index] = min(result);
gbest_pos = pos(g_index,:); 
gbest_res(1:FEs) = gbest_val;

obj_func_slope = zeros(num_g,1);
fri_best = (1:num_g1)'*ones(1,dimension);

%% Updateding examplers for the CL subpopulation(first time)
for i = 1:num_g1                                    
    fri_best(i,:) = i*ones(1,dimension);
    friend1 = ceil(num_g*rand(1,dimension));
    friend2 = ceil(num_g*rand(1,dimension));
    friend = (pbest_val(friend1) < pbest_val(friend2)).*friend1 + (pbest_val(friend1) >= pbest_val(friend2)).*friend2;
    toss = ceil(rand(1,dimension) - Pc(:,i)');
    if toss == ones(1,dimension)
       temp_index = randperm(dimension);           
       toss(1,temp_index(1)) = 0;
       clear temp_index;
    end
    fri_best(i,:) = (1-toss).*friend + toss.*fri_best(i,:);
    for d = 1:dimension
        fri_best_pos(i,d) = pbest_pos(fri_best(i,d),d);
    end
end

%% Updateding examplers for the DMS population
for i = 1:group_num                                
    group_id(i,:) = [((i-1)*group_ps + num_g1+1):i*group_ps + num_g1];
    pos_group(group_id(i,:)) = i;
    [gbestval(i),gbestid] = min(pbest_val(group_id(i,:)));    
    gbest(i,:) = pbest_pos(group_id(i,gbestid),:);            
end

n = 0;    mm = 0;

while n <= max_iteration && FEs<=max_FEs
    p = gbest_val;
    n = n + 1;
    Average_g = mean(result);                                       
    
    for i = 1:group_num
        Average_n(group_id(i,:)) = mean(result(group_id(i,:))); 
    end
    
    %% The CL subpopulation update velocity and position
    delta_g1 = (c(n).*rand(num_g1,dimension).*(fri_best_pos(1:num_g1,:) - pos(1:num_g1,:)));
    vel_g1 = Weight(n)*vel(1:num_g1,:) + delta_g1;
    vel_g1 = ((vel_g1 < v_min(1:num_g1,:)).*v_min(1:num_g1,:)) + ((vel_g1 > v_max(1:num_g1,:)).*v_max(1:num_g1,:)) + (((vel_g1 < v_max(1:num_g1,:))&(vel_g1 > v_min(1:num_g1,:))).*vel_g1);
    pos_g1 = pos(1:num_g1,:) + vel_g1;
    
    %% The DMS subpopulation update velocity and position
    for i = num_g1 + 1: num_g
        if Average_n(i) >= Average_g
            wx(i) = Weight1(n) + C;
            if wx(i)>0.99
                wx(i) = 0.99;
            end
        else
            wx(i) = Weight1(n) - C;
            if wx(i) < 0.20
                wx(i) = 0.20;
            end
        end
        delta_g2(i,:) = c1(n).*rand(1,dimension).*(pbest_pos(i,:) - pos(i,:)) + c2(n).*rand(1,dimension).*(gbest(pos_group(i),:) - pos(i,:));
        vel_g2(i,:) = wx(i)*vel(i,:) + delta_g2(i,:);
        
        vel_g2(i,:) = ((vel_g2(i,:) < v_min(i,:)).*v_min(i,:)) + ((vel_g2(i,:) > v_max(i,:)).*v_max(i,:)) + (((vel_g2(i,:) < v_max(i,:))&(vel_g2(i,:) > v_min(i,:))).*vel_g2(i,:));
        pos_g2(i,:) = pos(i,:) + vel_g2(i,:);
        keep_d = rand(1,dimension) < 0.5;
        pos_g2(i,:) = keep_d.*pbest_pos(i,:) + (1 - keep_d).*pos_g2(i,:); 
    end
    
    % The whole population
    pos_g2(1:num_g1,:) = [];
    vel_g2(1:num_g1,:) = [];
    pos = [pos_g1;  pos_g2];         
    vel = [vel_g1;  vel_g2];
    
    % Evaluate fitness
    for i=1:num_g
        if (sum(pos(i,:) > Xmax(i,:)) + sum(pos(i,:) < Xmin(i,:)) == 0)  
            result(i) = CostFunction(pos(i,:)',func_num); 
            FEs = FEs + 1;
            if  result(i) < pbest_val(i)
                pbest_pos(i,:) = pos(i,:);
                pbest_val(i) = result(i);
                obj_func_slope(i) = 0;               
            else
                obj_func_slope(i) = obj_func_slope(i) + 1;
            end
            
            if pbest_val(i) < gbest_val
                gbest_pos = pbest_pos(i,:);
                gbest_val = pbest_val(i);   
            end     
            gbest_res(FEs) = gbest_val; 
            if FEs >= max_FEs
                break;
            end
        end
    end
    
    %% update local gbest value and postion
    for i = num_g1+1:num_g
        if pbest_val(i) < gbestval(pos_group(i))
            gbest(pos_group(i),:) = pbest_pos(i,:);
            gbestval(pos_group(i)) = pbest_val(i);
        end
    end
    
    if gbest_val < p          
        flag = 0;
    else
        flag = flag + 1;
    end
    
    if flag >= gap1
        pt = gbest_pos;
        d1 = unidrnd(dimension);    randdata = 2 * rand(1,1) - 1;
        pt(d1) = pt(d1) + sign(randdata)*(UB(1) - LB(1))*normrnd(0,sig^2);     
        pt(find(pt(:)>UB(1))) = UB(1) * rand;                               
        pt(find(pt(:)<LB(1))) = LB(1) * rand;    
        cv = CostFunction(pt',func_num);  
        FEs = FEs+1;                                                                     
        if cv < gbest_val
            gbest_pos = pt;
            gbest_val = cv;
            flag=0;
        end
        gbest_res(FEs) = gbest_val;
        if FEs>=max_FEs
            break;
        end
    end
    sig = sigmax - (sigmax-sigmin)*(FEs/max_FEs);
    
    %% updating exemplar for the CL subpopulation
    for i = 1:num_g1
        if obj_func_slope(i) > 5
            fri_best(i,:) = i*ones(1,dimension);         
            friend1 = ceil(num_g1*rand(1,dimension));    
            friend2 = ceil(num_g1*rand(1,dimension));
            friend = (pbest_val(friend1) < pbest_val(friend2)).*friend1 + (pbest_val(friend1) >= pbest_val(friend2)).*friend2;
            toss = ceil(rand(1,dimension) - Pc(:,i)');
            
            if toss == ones(1,dimension)
                temp_index = randperm(dimension);
                toss(1,temp_index(1)) = 0;
                clear temp_index;
            end
            
            fri_best(i,:) = (1 - toss).*friend + toss.*fri_best(i,:);
            for d = 1:dimension
                fri_best_pos(i,d) = pbest_pos(fri_best(i,d),d);
            end
            obj_func_slope(i) = 0;
        end
    end 
    
    %% Regrouping the sub-swarm particles for the DMS subpopulation
    if mod(n,100)==0
        options = optimset('LargeScale','off','MaxFunEvals',L_FEs,'Display','off');
        [tmp,tmpid]=sort(gbestval);
        for k=1:L_num
            [x,fval,exitflag,output] = fminunc(fhd1,gbest(tmpid(k),:)',options,func_num);             
            FEs=FEs+L_FEs;
            if fval<gbestval(tmpid(k))
                [gbestval(tmpid(k)),gbestid] = min(pbest_val(group_id(tmpid(k),:)));
                pbest_pos(group_id(tmpid(k),gbestid),:) = x'; 
                pbest_val(group_id(tmpid(k),gbestid)) = fval;
                gbest(tmpid(k),:) = x'; 
                gbestval(tmpid(k)) = fval;
            end
            gbest_res(1,FEs-L_FEs+1:FEs) = gbest_val;                    
            [gbest_val,index2] = min([gbestval gbest_val]);                
            gbest_postmp = [gbest;gbest_pos];
            gbest_pos = gbest_postmp(index2,:);

            if FEs >= max_FEs
                break;
            end
        end
    end
    
    if mod(n,R(n)) == 0                     
        rc = randperm(num_g2) + num_g1;
        group_id = []; gbest = []; gbestval = [];
        for i = 1:group_num
            group_id(i,:) = rc(((i-1)*group_ps + 1):i*group_ps);
            pos_group(group_id(i,:)) = i;
            [gbestval(i),gbestid] = min(pbest_val(group_id(i,:)));
            gbest(i,:) = pbest_pos(group_id(i,gbestid),:);
        end
        
        for i=num_g1 + 1: num_g  
            posa(i,:)=pos(i,:);
        end
        for i = 1: num_g2 
            pos(i+num_g1,:)=posa(rc(i),:);
        end
       
        for i=num_g1 + 1: num_g  
            vela(i,:)=vel(i,:);
        end
        for i = 1: num_g2 
            vel(i+num_g1,:)=vela(rc(i),:);
        end       
    end 
    
    if FEs>=max_FEs
        break;
    end
    if (n==max_iteration) && (FEs<max_FEs)
        n=n-1;
    end
     
end  % end while

FEs = max_FEs; 
Best.val = gbest_val;
Best.pos = gbest_pos;

end    
       
  
