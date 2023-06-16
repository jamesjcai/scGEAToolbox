function [X,genelist]=r_SoupX(selpath)
%Run SoupX decontamination
%
% https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
%
% see also: run.r_decontX

oldpth=pwd();
X=[];
genelist=[];

[isok,msg]=commoncheck_R('R_SoupX');
if ~isok, error(msg); return; end
if nargin<1
    selpath = uigetdir([],'CellRanger outs Folder');
end
if ~(ischar(selpath) && exist(selpath,'dir')==7)
    return;
end
if exist('./input.txt','file'), delete('./input.txt'); end
if exist('./output.h5','file'), delete('./output.h5'); end

writematrix(selpath,'input.txt');
Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);
if exist('./output.h5','file')
    X=h5read('output.h5','/X');
    try 
        disp('Reading genes...')
        genelist=i_getgenelist(selpath);
    catch
        cd(oldpth);
        return;
    end
end
if exist('./input.txt','file'), delete('./input.txt'); end
if exist('./output.h5','file'), delete('./output.h5'); end
cd(oldpth);
end




function [genelist]=i_getgenelist(inpath)
    selpath=fullfile(inpath,'filtered_feature_bc_matrix');    
    ftdone=false;
    if ~ftdone
        ftfname=fullfile(selpath,'features.tsv');
        zftfname=fullfile(selpath,'features.tsv.gz');
        if ~exist(ftfname,'file')
            if ~exist(zftfname,'file')
                % error('No features.tsv file.');
            else
                [~,nametxt]=fileparts(zftfname);
                fprintf('Unzipping %s.gz...\n',nametxt);
                gunzip(zftfname);
                ftdone=true;
            end
        else
            ftdone=true;
        end
    end

    if ~ftdone
        ftfname=fullfile(selpath,'genes.tsv');
        zftfname=fullfile(selpath,'genes.tsv.gz');
        if ~exist(ftfname,'file')
            if ~exist(zftfname,'file')
                % error('No features.tsv file.');
            else
                [~,nametxt]=fileparts(zftfname);
                fprintf('Unzipping %s.gz...\n',nametxt);
                gunzip(zftfname);
                ftdone=true;
            end
        else
            ftdone=true;
        end
    end
    if ftdone
        T=readtable(ftfname,'ReadVariableNames',false,...
            'filetype','text','Delimiter',{'\t',',',' ',';','|'});
        genelist=string(T.Var2);
    end
end
