    function obj = estimatepotency(obj,speciesid,forced)
        if nargin<3, forced=false; end
        if nargin<2, speciesid=[]; end
        if forced || sum(strcmp('cell_potency',obj.list_cell_attributes))==0
            if isempty(speciesid)               
                % speciesid=input('Species: 1=human,2=mouse >>');
                answer = questdlg('Which species?','Select Species','Mouse','Human','Mouse');

                if strcmp(answer,'Human')
                    speciesid=1;
                elseif strcmp(answer,'Mouse')
                    speciesid=2;
                else
                    return;
                end
            elseif ischar(speciesid)
                [y,idx]=ismember(lower(speciesid),{'human','mouse'});
                if y, speciesid=idx; end
            elseif isstring(speciesid)
                [y,idx]=ismember(lower(speciesid),["human","mouse"]);
                if y, speciesid=idx; end
            elseif isnumeric(speciesid)
                if ~ismember(speciesid,[1 2])
                    error('SPECIESID should be 1 (human) or 2 (mouse)');
                end
            end
            r=sc_potency(obj.X,obj.g,speciesid);
            obj.list_cell_attributes=[obj.list_cell_attributes,...
                {'cell_potency',r}];
            disp('cell_potency added.');
        else
            disp('cell_potency existed.');
        end
    end
