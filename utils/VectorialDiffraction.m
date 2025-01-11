function [Ex,Ey,Ez] = VectorialDiffraction(Obj,Beam,Scope,FocusingMethod,number_u,number_v)
%VECTORIALDIFFRACTION calculate vectorial diffraction
%
% *************************************************************************
% originates from Dr. Hao
% optimized by LIU Xin
% liuxin24@hku.hk
% Apr.24, 2021

switch FocusingMethod
    case 'int_ap'  % GT
        [Ex,Ey,Ez] = VectorialDiffraction_GT_AP(Obj,Beam,Scope);
    case 'smm_ap'  % our method
        [Ex,Ey,Ez] = VectorialDiffraction_ASR_AP(Obj,Beam,Scope,number_u,number_v);
    case 'czt_ap'  % Ctrl
        [Ex,Ey,Ez] = VectorialDiffraction_CZT_AP(Obj,Beam,Scope);
    otherwise
        error('No such a mode, please try again!');
end

end