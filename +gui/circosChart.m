classdef circosChart
% @author : slandarer
% 公众号  : slandarer随笔 
% 知乎    : hikari

    properties
        ax,arginList={'ColorOrder','ClassName','PartName'}
        ColorOrder=[80,118,169;226,144,50;127,167,58;242,86,54;126,109,167;
                    196,98,37;74,148,189;255,182,46;161,86,144;134,138,33;
                    240,73,53;90,123,207;254,147,44;186,79,115;35,170,102]./255;
        ClassName,PartName
        Data,Class,indexInClass,colorSet={[]}
        classSet,classNum,classSize,classRatio,classTheta
        lineHdl,partLabelHdl,classLabelHdl,scatterHdl
    end

    methods
        function obj=circosChart(Data,Class,varargin)
            obj.Data=Data;
            obj.Class=Class(:);
            obj.classSet=unique(Class);
            obj.classNum=length(obj.classSet);

            obj.indexInClass=zeros(length(obj.Class),1);
            % 计算比例
            for i=1:obj.classNum
                tClassBool=obj.classSet(i)==obj.Class;
                tCumsumBool=cumsum(tClassBool);
                obj.classSize(i)=sum(tClassBool);
                obj.indexInClass(tClassBool)=tCumsumBool(tClassBool);
            end
            obj.classRatio=obj.classSize./sum(obj.classSize);
            disp(char([64 97 117 116 104 111 114 32 58 32 115 108 97 110 100 97 114 101 114]))
            obj.ColorOrder=[obj.ColorOrder;rand([obj.classNum,3])];
            for i=1:size(obj.Data,1)
                obj.PartName{i}='';
            end
            for i=1:obj.classNum
                obj.ClassName{i}=['Class ',num2str(i)];
            end
            % 获取其他数据
            for i=1:2:(length(varargin)-1)
                tid=ismember(obj.arginList,varargin{i});
                if any(tid)
                obj.(obj.arginList{tid})=varargin{i+1};
                end
            end
        end

        function obj=draw(obj)
            obj.ax=gca;hold on
            obj.ax.XLim=[-1.2,1.2];
            obj.ax.YLim=[-1.2,1.2];
            obj.ax.XTick=[];
            obj.ax.YTick=[];
            obj.ax.XColor='none';
            obj.ax.YColor='none';
            obj.ax.PlotBoxAspectRatio=[1,1,1];

            % 调整初始界面大小
            fig=obj.ax.Parent;
            % if max(fig.Position(3:4))<600
            %     fig.Position(3:4)=1.8.*fig.Position(3:4);
            %     fig.Position(1:2)=fig.Position(1:2)./3;
            % end

            sepTheta=2/30/obj.classNum;
            cumTheta=[0,28/30*cumsum(obj.classRatio)];

            if isempty(obj.PartName{1})&&isempty(obj.PartName{2})
                tdis=1.1;
            else
                tdis=1.22;
            end
            % 计算每一类中每一个元素的角度
            for i=1:obj.classNum
                obj.classTheta(i).T=linspace(sepTheta*i+cumTheta(i),sepTheta*i+cumTheta(i+1),obj.classSize(i)).*2.*pi;
                obj.scatterHdl(i)=scatter(cos(obj.classTheta(i).T),sin(obj.classTheta(i).T),30,'filled','CData',obj.ColorOrder(i,:),'MarkerEdgeColor','none');
                
                CTi=mean(obj.classTheta(i).T);
                rotation=CTi/pi*180;
                if rotation>0&&rotation<180
                    obj.classLabelHdl(i)=text(cos(CTi).*tdis,sin(CTi).*tdis,obj.ClassName{i},'FontSize',14,'FontName','Arial',...
                    'HorizontalAlignment','center','Rotation',-(.5*pi-CTi)./pi.*180);
                else
                    obj.classLabelHdl(i)=text(cos(CTi).*tdis,sin(CTi).*tdis,obj.ClassName{i},'FontSize',14,...
                    'HorizontalAlignment','center','Rotation',-(1.5*pi-CTi)./pi.*180);
                end
            end
            % 绘制文字
            for i=1:size(obj.Data,1)

                Ci=obj.Class(i);Pi=obj.indexInClass(i);
                Ti=obj.classTheta(Ci).T(Pi);
                rotation=Ti/pi*180;
                 
                if rotation>90&&rotation<270
                    rotation=rotation+180;
                    obj.partLabelHdl(i)=text(cos(Ti).*1.03,sin(Ti).*1.03,obj.PartName{i},'Rotation',rotation,'HorizontalAlignment','right','FontSize',8);
                else
                    obj.partLabelHdl(i)=text(cos(Ti).*1.03,sin(Ti).*1.03,obj.PartName{i},'Rotation',rotation,'FontSize',8);
                end
            end

            % 计算类与类之间的渐变色
            t2=linspace(0,1,200);t1=1-t2;
            for i=1:obj.classNum
                for j=1:obj.classNum
                    C1=obj.ColorOrder(i,:);
                    C2=obj.ColorOrder(j,:);
                    obj.colorSet{i,j}=uint8([t1.*C1(1)+t2.*C2(1);
                        t1.*C1(2)+t2.*C2(2);
                        t1.*C1(3)+t2.*C2(3)
                        ones(1,200).*.6].*255);
                end
            end
            % 画线并赋予颜色
            for i=1:size(obj.Data,1)
                for j=1:(i-1)
                    if obj.Data(i,j)>0
                        Ci=obj.Class(i);Pi=obj.indexInClass(i);
                        Cj=obj.Class(j);Pj=obj.indexInClass(j);
                        Ti=obj.classTheta(Ci).T(Pi);
                        Tj=obj.classTheta(Cj).T(Pj);
                        Xij=[cos(Ti),0,cos(Tj)]';
                        Yij=[sin(Ti),0,sin(Tj)]';
                        XYb=bezierCurve([Xij,Yij],200);
                        obj.lineHdl(i,j)=plot(XYb(:,1),XYb(:,2),'-','LineWidth',1);
                    end
                end
            end
            pause(1e-16)
            for i=1:size(obj.Data,1)
                for j=1:(i-1)
                    if obj.Data(i,j)>0 
                        Ci=obj.Class(i);
                        Cj=obj.Class(j);
                        set(get(obj.lineHdl(i,j),'Edge'),'ColorBinding','interpolated','ColorData',obj.colorSet{Ci,Cj})
                    end
                end
            end
            % 贝塞尔函数
            function pnts=bezierCurve(pnts,N)
                t=linspace(0,1,N);
                p=size(pnts,1)-1;
                coe1=factorial(p)./factorial(0:p)./factorial(p:-1:0);
                coe2=((t).^((0:p)')).*((1-t).^((p:-1:0)'));
                pnts=(pnts'*(coe1'.*coe2))';
            end
        end
        % 设置线除了颜色的其他属性
        function setLine(obj,varargin)
             for i=1:size(obj.Data,1)
                for j=1:(i-1)
                    if obj.Data(i,j)>0
                        set(obj.lineHdl(i,j),varargin{:})
                    end
                end
             end
        end
        % 设置线颜色
        function setColor(obj,N,color)
            obj.ColorOrder(N,:)=color;
            t2=linspace(0,1,200);t1=1-t2;
            for i=1:obj.classNum
                set(obj.scatterHdl(i),'CData',obj.ColorOrder(i,:))
                for j=1:obj.classNum
                    C1=obj.ColorOrder(i,:);
                    C2=obj.ColorOrder(j,:);
                    obj.colorSet{i,j}=uint8([t1.*C1(1)+t2.*C2(1);
                        t1.*C1(2)+t2.*C2(2);
                        t1.*C1(3)+t2.*C2(3)
                        ones(1,200).*.6].*255);
                end
            end
            for i=1:size(obj.Data,1)
                for j=1:(i-1)
                    if obj.Data(i,j)>0 
                        Ci=obj.Class(i);
                        Cj=obj.Class(j);
                        set(get(obj.lineHdl(i,j),'Edge'),'ColorBinding','interpolated','ColorData',obj.colorSet{Ci,Cj})
                    end
                end
            end
        end

        % 设置标签
        function setPartLabel(obj,varargin)
            for i=1:size(obj.Data,1)
                set(obj.partLabelHdl(i),varargin{:});
            end
        end

        function setClassLabel(obj,varargin)
            for i=1:obj.classNum
                set(obj.classLabelHdl(i),varargin{:});
            end
        end
    end
% @author : slandarer
% 公众号  : slandarer随笔 
% 知乎    : hikari
end



%{
% demo 1 随机数据
rng(1)


% 生成随机200x200对称0-1矩阵
Data=rand(200,200)>.992;
sum(sum(Data))
Data=(Data+Data')>0;
% 生成200x1随机分类编号
Class=randi([1,5],[200,1]);

for i=1:200
    partName{i}=[num2str(Class(i)),'-',num2str(i)];
end
className={'AAAAA','BBBBB','CCCCC','DDDDD','EEEEE'};

% CC=circosChart(Data,Class);
CC=circosChart(Data,Class,'PartName',partName,'ClassName',className);
CC=CC.draw();

CC.setPartLabel('Color',[0,0,.8],'FontName','Cambria')
CC.setClassLabel('Color',[.8,0,0],'FontName','Cambria','FontSize',25)
%}