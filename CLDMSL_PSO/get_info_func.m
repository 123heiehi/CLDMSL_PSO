function [dimension , Search_Range] = get_info_func(what_func)

    switch what_func
        case 1
            dimension   = 9;
            Search_Range       	= [-8192, 8192];
        case 2
            dimension   = 16;
            Search_Range       	= [-16384, 16384];
        case 3
            dimension   = 18;
            Search_Range       	= [-4,4];
        case 4
            dimension   = 10;
            Search_Range       	= [-100,100];
        case 5
            dimension   = 10;
            Search_Range       	= [-100,100];
        case 6
            dimension   = 10;
            Search_Range       	= [-100,100];
        case 7
            dimension   = 10;
            Search_Range       	= [-100,100];
        case 8
            dimension   = 10;
            Search_Range       	= [-100,100];
        case 9
            dimension   = 10;
            Search_Range       	= [-100,100];
        case 10
            dimension   = 10;
            Search_Range       	= [-100,100];
    end
end