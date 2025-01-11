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
% updated by Xin Liu
% email: liuxin24@hku.hk
% Mar.19, 2020

% Previous percentage number.
persistent i_prev_prct;

% Current percentage number.
i_prct = i ./ n * 100;

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
    S_finish = ' Done!\n';
    fprintf(S_finish);
    clear i_prev_prct;
%     fprintf(getBackspaceStr(length(msg)+2+length(S)+length(S_finish)));
end

end

function S = getPrctStr(prct)
% control the precision of the percentage
S = sprintf('%.2f%% %s',prct,getDotStr(prct));

if prct < 10
    S = ['  ',S];
elseif prct < 100
    S = [' ',S];
end
end

function S = getDotStr(prct)
S = repmat(' ',1,10);
S(1:floor(prct/10)) = '>';
S = ['[',S,'].'];
end

function S = getBackspaceStr(N)
S = repmat('\b',1,N);
end