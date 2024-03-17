classdef TestMDS< handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    methods(Static)
        function cityMap=NorthAmericanCities
            cities = {'Edmonton', 'Winnipeg', 'Portland', ...
                'Atlanta','Chicago','Denver',...
                'Houston','L.A.','Miami',...
                'New York','SF bay area','Seattle','Washington DC', ...
                'Vancouver', 'Toronto', 'Montreal', 'St. John''s', ...
                'Dallas', 'Halifax', 'Juneau',...
                'Minneapolis', 'Quebec city', 'Memphis',...
                'Tampa Bay', 'Mexico city'};
            populations=[1321426 778489 2389228 ...%from Google
                4789700 9512999 2814330 ...
                6313158 10176000 5564635 ...
                19800000   7000000 3733580 6131977 ...
                2463431 5928040 4098927 205955 ...
                6426214 403390 100000 ...
                3551036 800296 1341746 ...
                2823938 21200000];
            latitudeLongtitude=...%from Google
                [53.5444 113.4909; 49.8951 97.1384; 45.5231 122.6765;...
                33.7490 84.3880; 41.8781 87.6298; 39.7392 104.9903;...
                29.7604 95.3698; 34.0522 118.2437; 25.7617 80.1918;...
                40.7128 74.0060; 37.7749 122.4194;47.6062 122.3321; 38.9072 77.0369;...
                49.2827 123.1207;43.6532 79.3832; 45.5017 73.5673; 47.5615 52.7126;...
                32.7767 96.7970;44.6488 63.5752; 58.3019 134.4197;...
                44.9778 93.2650; 46.8139 71.2080;35.1495 90.0490; ...
                27.7634 82.5437; 19.4326 99.1332;];

            cityMap=MDS(cities, latitudeLongtitude, populations);
            [fig, ~]=cityMap.newFig('North American cities', ...
                false, {gcf, 'east', true}, 'TestMDS.fig');
            p=get(fig, 'position');
            set(fig, 'position', [p(1)*.8 p(2) p(3)*1.4 p(4)*1.4]);
            ax=get(fig, 'currentaxes');
            if isempty(ax)
                ax=gca;
            end
            %orient map west-east and north-south
            cityMap.XY=TestMDS.flip(cityMap.XY, 1);
            cityMap.XY=TestMDS.flip(cityMap.XY, 2);
            cityMap.setNames('Coordinate', 'lattitude/longtitude', 'North American Cities');
            cityMap.setMeasurementNames({'Lat','Long'});
            cityMap.setSeeTextProperty('TestMDS');
            %cityMap.setSymbolDimension('Lat');
            sp=zeros(1, cityMap.nThings);
            sp(1:3:end)=3;
            sp(2:3:end)=2;
            %cityMap.setSeePositions(sp);
            cityMap.see;
            cityMap.see(2);
            cityMap.see(3);
            cityMap.setMouseEar(@notifyMouse);

            function notifyMouse(~, isMotion, idx,x, y)
                if ~isMotion
                    msg(['Clicked on ' cities{idx}]);
                end
            end
        end

        function m=flip(m, d)
            mx=max(m(:,d));
            mn=min(m(:,d));
            m(:,d)=mx-(m(:,d)-mn);
        end


    end
end