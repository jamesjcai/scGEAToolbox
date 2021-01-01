    function obj = estimatepotency(obj)
        if sum(strcmp('cell_potency',obj.list_cell_attributes))==0;
            idx=input('Species: 1=human,2=mouse >>');
            r=sc_potency(obj.X,obj.g,idx);
            obj.list_cell_attributes=[obj.list_cell_attributes,...
                {'cell_potency',r}];
            disp('cell_potency added.');
        else
            disp('cell_potency existed.');
        end
    end
