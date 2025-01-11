function Eout = ScalarDiffraction(Beam,Scope,DiffractionMethod,number_u,number_v)
%SCALARDIFFRACTION calculate scalar diffraction
%
% *************************************************************************
% LIU Xin
% liuxin24@hku.hk
% Apr.23, 2021

switch DiffractionMethod
    case 'asap'  % GT
        Eout = ScalarDiffraction_GT_AP(Beam,Scope);
    case 'assmm'  % our method
        Eout = ScalarDiffraction_ASR_AP(Beam,Scope,number_u,number_v);
    case 'asczt'  % Ctrl
        Eout = ScalarDiffraction_CZT_AP(Beam,Scope);
    otherwise
        error('No such a mode, please try again!');
end

end