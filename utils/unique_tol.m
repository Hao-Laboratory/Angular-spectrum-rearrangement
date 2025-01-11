function [C,IC] = unique_tol(A,num)
%UNIQUE_TOL obtain the spatial frequency and the corresponding index after
%angular spectrum merging without weight calculation
%
% INPUT********************************************************************
% A: original spatial frequency
% num: number of target samples in the frequency domain after merging
% OUTPUT*******************************************************************
% C: spatial frequency after merging without weight calculation
% IC: the corresponding index
% 
% HU Yiwen
% huyw@zju.edu.cn

[A, ~, idx_A] = unique(A);
N = length(A(:));

A_temp2 = A;
ref = nan;
C_A_temp2 = [];

if N>num
    while (length(A_temp2)>num)
        idxu = length(A_temp2) - num;
        dA = diff(A_temp2);
        dA2 = sort(dA);
        if dA2(idxu)<=dA2(idxu)/max(abs(A_temp2(:)))*max(abs(A_temp2(:)))
            A_temp1 = A_temp2;  % record last frequency
            A_temp2 = uniquetol(A_temp2,dA2(idxu)/max(abs(A_temp2(:))));
        else  % avoid numerical error
            for a = sort(A_temp2(:)).'
                if ~(a-ref <= dA2(idxu))
                    ref = a;
                    C_A_temp2(end+1) = ref;
                end
            end
                A_temp1 = A_temp2;
                A_temp2 = C_A_temp2.';
                C_A_temp2 = [];
                ref = nan;
        end
    end

   % avoid excessive merging
   if length(A_temp2) < num
        idxuu = find(dA == dA2(idxu));
        [~, idx_A_temp2] = ismember(A_temp2,A_temp1);
        idx_A_temp2_pre = idx_A_temp2(ismember(idx_A_temp2,idxuu));
        M = num - length(A_temp2);
        idx_A_temp2_add = idx_A_temp2_pre + 1;
        A2_add = A_temp1(idx_A_temp2_add(1:M));
        A_temp2 = [A_temp2;A2_add];
        A_temp2 = sort(A_temp2);
   end

   C = A_temp2;

   idx_eff = 1:num;
   [~, idx_C] = ismember(C,A);
   idx_C = reshape(idx_C,[1,length(idx_C)]);
   repNum = diff([idx_C,N+1]);

   ic = repelem(idx_eff,repNum);
   IC = ic(idx_A);
else
    C = A;
    ic = 1:N;
    IC = ic(idx_A);
end