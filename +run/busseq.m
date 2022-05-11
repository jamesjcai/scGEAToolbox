function [Y0,Y1]=busseq(X0,X1)
% BUSseq (Batch Effects Correction with Unknown Subtypes for scRNA seq)

% BUSseq is an R tool to correct batch effects in single-cell RNA-seq data. 
% The BUSseq algorithm can also cluster cell types, correct for 
% overdispersion, cell-specific sequencing depth, and dropout events.
%
% NCELLS=2000; NGENES=400;
% X0=nbinrnd(20,0.98,NGENES,NCELLS);
% X0=X0(:,sum(X0)>120);
% X1=nbinrnd(20,0.98,NGENES,NCELLS);
% X1=X1(:,sum(X1)>120);


if nargin<2, error('USAGE: [Y0,Y1]=run.busseq(X0,X1)'); end
if nargout<2, error('USAGE: [Y0,Y1]=run.busseq(X0,X1)'); end
Y0=[]; Y1=[];

isdebug=true;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_BUSseq');
if ~isok, error(msg); end



tmpfilelist={'input.h5','output.h5'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

h5create('input.h5', '/X0', size(X0));
h5write('input.h5', '/X0', X0);
h5create('input.h5', '/X1', size(X1));
h5write('input.h5', '/X1', X1);

writematrix(X0,'input0.csv');
writematrix(X1,'input1.csv');
% for i = 1:50
%     textwaitbar(i, 100, 'This may take a few minutes. Please wait');
%     pause(0.01);
% end
pkg.RunRcode('Copy_of_script.R');
% for i = 51:100
%     textwaitbar(i, 100, 'This may take a few minutes. Please wait');
%     pause(0.01);
% end
if exist('output.h5','file')
    %Y0=readmatrix('output0.csv');
    %Y1=readmatrix('output1.csv');
    Y0=h5read("output.h5",'/Y0');
    Y1=h5read("output.h5",'/Y1');
end
% if exist('input0.csv','file'), delete('input0.csv'); end
% if exist('output0.csv','file'), delete('output0.csv'); end
% if exist('input1.csv','file'), delete('input1.csv'); end
% if exist('output1.csv','file'), delete('output1.csv'); end
tmpfilelist={'input.h5','output.h5'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

cd(oldpth);

end


function textwaitbar(i, n, msg)
% A command line version of waitbar.
% Usage:
%   textwaitbar(i, n, msg)
% Input:
%   i   :   i-th iteration.
%   n   :   total iterations.
%   msg :   text message to print.
%
% Date      : 05/23/2019
% Author    : Xiaoxuan He   <hexxx937@umn.edu>
% Institute : University of Minnesota
%
% Previous percentage number.
persistent i_prev_prct;
% Current percentage number.
i_prct = floor(i ./ n * 100);
% Print message when counting starts.
if isempty(i_prev_prct) || i_prct < i_prev_prct
    i_prev_prct = 0;
    S_prev = getPrctStr(i_prev_prct);
    
    fprintf('%s: %s',msg, S_prev);
end
% Print updated percentage.
if i_prct ~= i_prev_prct
    S_prev = getPrctStr(i_prev_prct);
    fprintf(getBackspaceStr(numel(S_prev)));
    
    S = getPrctStr(i_prct);
    fprintf('%s', S);
    
    i_prev_prct = i_prct;
end
% Clear percentage variable.
if i_prct == 100
    fprintf(' Done.\n');
    clear i_prev_prct;
end
end
function S = getPrctStr(prct)
S = sprintf('%d%%  %s',prct,getDotStr(prct));
if prct < 10
    S = ['  ',S];
elseif prct < 100
    S = [' ',S];
end
end
function S = getDotStr(prct)
S = repmat(' ',1,10);
S(1:floor(prct/10)) = '.';
S = ['[',S,']'];
end
function S = getBackspaceStr(N)
S = repmat('\b',1,N);
end

