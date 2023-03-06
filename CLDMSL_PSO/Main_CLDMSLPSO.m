%----------- CLDMSL_PSO for the  the 100-Digit Challenge --------------------
%---------------------- Written by Rui Wang---------------------

clear all ;  format longG;
method_name  = 'CLDMSL_PSO';
disp(['Hello! Please wait . . . ' method_name]);
diary (['DiaryFile '    method_name]);
CostFunction = @(pop,the_func)     cec19_func(pop,the_func);
list_func    = [1:10];
numb_func    = length(list_func);
for i        = 1 : numb_func
    tic
    the_func = list_func(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       U S E R    D E F I N I T I O N
            repeated        = 30;
    %--------------- Initial Parameters of CLDMSL_PSO---------------------------
            para.PopSize    = 50;

    %----------------------------------------------------------------------
            [dimension,Search_Range] = get_info_func(the_func);
            Info.f_star         = 1.000000000;
            Info.func_num       = the_func;
            Info.max_FEs        = 300000;   % 1e9;
            Info.dimension      = dimension;
            Info.Search_Range   = Search_Range;
    %       E N D    O F     U S E R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor rep  = 1 : repeated
        [Best, FEs]     = CLDMSL_PSO(Info, para, CostFunction);
%         [Best , array_digit , FEs , Mig]     = SOMA_T3A(Info,SOMApara,CostFunction);
        % ----------- Record Temp-Values ----------------------------------
        error                   = Best.val - Info.f_star;
        path_error(rep)         = error;

    end     % END REPEATED  (parfor Loop)
    %----------- Record Data ----------------------------------------------
    data(the_func).time             = toc;
    data(the_func).path_error       = path_error;
    data(the_func).mean_error       = mean(data(the_func).path_error);
    data(the_func).std_error        = std(data(the_func).path_error);

% end    %  E N D   F U N C
disp('===================================================================');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(method_name,'data');
end    %  E N D   F U N C
diary off
%------------------------ END OF EVERYTHING -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%