clear;

%%Load data
A = table2array(load_data_from_spectrometer_file('Data/Real_images/Calculated ASCII RAD'));

%% Take the mean of the spectra
A=A';
Z=[];
for i=1:10:size(A,1)-9
    tmp=A(i:i+9,:);
    Z=[Z; mean(tmp)];
end

%% Plotting the spectra

figure;
axis([50 650 0 1]);
plot(Z(:,50:650)');

color=Z(2,50:650)./Z(1,50:650);

