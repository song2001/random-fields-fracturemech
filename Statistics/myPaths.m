% Vincente Pericoli
% UC Davis
% random-fields-fracturemech
% 8 Dec 2015
%

function dir_str = myPaths(lib_str)
%MYPATHS: Given a library name (string), returns the directory (string)
%         where it exists.
%
%Defined libraries:
%   VGPy_databases
%   fem-interp
%


switch lower(lib_str)
    case 'vgpy_databases'
        dir_str = 'C:\Temp\VGPy_Databases';
    case 'fem-interp'
        dir_str = 'C:\Users\Vince Pericoli\Documents\GitHub\fem-interp';
end


end


