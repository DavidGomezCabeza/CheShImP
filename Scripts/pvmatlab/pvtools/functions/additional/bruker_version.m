function bruker_version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2012
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

superclass=BrukerDataSuperclass;
str=sprintf('Package-Name: \t\t %s', superclass.PackageParameters.MatlabPackageName);
disp(str);
str=sprintf('Package-Version: \t %s', superclass.PackageParameters.MatlabPackageVersion);
disp(str);


end

