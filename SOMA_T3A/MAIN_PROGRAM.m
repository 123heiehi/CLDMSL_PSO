%----------- SOMA T3A for the  the 100-Digit Challenge --------------------
%-------- Written by Quoc Bao DIEP (diepquocbao@gmail.com) ----------------

clear all ;  format longG;
method_name  = 'SOMA_T3A';
disp(['Hello! Please wait . . . ' method_name]);
% diary (['DiaryFile '    method_name]);
CostFunction = @(pop,the_func)     cec19_func(pop,the_func);
list_func    = [1:10];
numb_func    = length(list_func);
for i        = 1 : numb_func
    tic
    the_func = list_func(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       U S E R    D E F I N I T I O N
            repeated            = 30;
    %--------------- Initial Parameters of SOMA ---------------------------
            SOMApara.PopSize    = 1500;
            SOMApara.N_jump     = 100;
			SOMApara.m     		= 5;
			SOMApara.n     		= 4;
			SOMApara.k     		= 5;
    %----------------------------------------------------------------------
            [dimension,Search_Range] = get_info_func(the_func);
            Info.f_star         = 1.000000000;
            Info.the_func       = the_func;
            Info.FEs_Max        = 300000;   % 1e9;
            Info.dimension      = dimension;
            Info.Search_Range   = Search_Range;
    %       E N D    O F     U S E R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor rep  = 1 : repeated
        [Best , array_digit , FEs , Mig]     = SOMA_T3A(Info,SOMApara,CostFunction);
        % ----------- Record Temp-Values ----------------------------------
        error                   = Best.Value - Info.f_star;
        path_error(rep)         = error;
        Mig_stop(rep)           = Mig;
        FEs_stop(rep)           = FEs.stop;
        FEs_reach_digit(:,rep)  =[FEs.reach_digit , FEs.stop]';
        matrix_digit(:,rep)     = array_digit';
        disp([FEs.reach_digit , FEs.stop , rep]');
    end     % END REPEATED  (parfor Loop)
    %----------- Record Data ----------------------------------------------
    data(the_func).time             = toc;
    data(the_func).path_error       = path_error;
    data(the_func).mean_error       = mean(data(the_func).path_error);
    data(the_func).std_error        = std(data(the_func).path_error);
    data(the_func).Mig_stop         = Mig_stop;
    data(the_func).FEs_stop         = FEs_stop;
    data(the_func).FEs_reach_digit  = FEs_reach_digit;
    data(the_func).matrix_digit     = matrix_digit;
    %----------- Display Value --------------------------------------------
    %clc
    disp('===================================================================');
    disp(['               ',method_name,'                            The Last']);
    disp('Function   Mean     (Std Dev)          Time(s)  Migrations     FEs');
    disp('-------------------------------------------------------------------');
%     for k   = 1:i
%         sub_func   = list_func(k);
%         mean_error = mean(data(sub_func).path_error);
%         std_error  =  std(data(sub_func).path_error);
% 
%         fprintf('   %2.0f      %8.2e (%8.2e)        %6.1f     %6.0f     %7.0f\n',...
%         sub_func,mean_error,std_error, data(sub_func).time, data(sub_func).Mig_stop(end), data(sub_func).FEs_stop(end));
%     end
    %----------- END Display Value ----------------------------------------
% end    %  E N D   F U N C
disp('===================================================================');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(method_name,'data');
data.FEs_reach_digit
end    %  E N D   F U N C
diary off
%------------------------ END OF EVERYTHING -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%