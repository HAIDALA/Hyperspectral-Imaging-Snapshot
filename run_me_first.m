
fprintf('###############################################################################\n');
fprintf('###                                                    ########################\n');
fprintf('###                Install Hyperspectral Imaging Toolbox first              ###\n');
fprintf('###                                                    ########################\n');
fprintf('###############################################################################\n');
addpath(pwd);

cd WNMFLibrary/;
addpath(genpath(pwd));
cd ..;


cd Data/
addpath(genpath(pwd));
cd ..;

cd Evaluate/
addpath(genpath(pwd));
cd ..;

cd Functions/
addpath(genpath(pwd));
cd ..;

cd Methods/
addpath(genpath(pwd));
cd ..;

cd Metrics/
addpath(genpath(pwd));
cd ..;

cd Results/
addpath(genpath(pwd));
cd ..;

cd Other_Methods/
addpath(genpath(pwd));
cd ..;
