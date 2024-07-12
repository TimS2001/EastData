clear
clc
close all

BC501A_Energy_calibration = [2.10803 46.08857];

currentPath = pwd;
disp(currentPath);
filename_BC501A_N5510 = fullfile(currentPath, 'exp_data', 'BC501A_N5510.mat');                 %需要处理的文件名
BC501A_N5510 = load(filename_BC501A_N5510);
BC501A_N5510_Time_Stamp_value = BC501A_N5510.max_Time_Stamp_value;
BC501A_N5510_spectrum_exp = BC501A_N5510.draw_Recoil_Proton;
figure('name','BC501A 5.51MeV Neutron recoil proton spectrum');
plot( BC501A_N5510_spectrum_exp(:,1)*BC501A_Energy_calibration(1)+BC501A_Energy_calibration(2)/1000, BC501A_N5510_spectrum_exp(:,2),'color','red','Linewidth',2);
title('BC501A 5.51MeV Neutron recoil proton spectrum','FontSize',20,'Color','black','FontWeight','bold');
xlabel('Equivalent electronic energy(MeVee)','FontSize',20,'Color','black','FontWeight','bold');
ylabel('Counts','FontSize',20, 'Color','black','FontWeight','bold');


EJ301_Energy_calibration  = [2.17838 22.79964];

currentPath = pwd;
disp(currentPath);
filename_EJ301_N5510 = fullfile(currentPath, 'exp_data', 'EJ301_N5510.mat');                 %需要处理的文件名
EJ301_N5510 = load(filename_EJ301_N5510);
EJ301_N5510_Time_Stamp_value = EJ301_N5510.max_Time_Stamp_value;
EJ301_N5510_spectrum_exp = EJ301_N5510.draw_Recoil_Proton;
figure('name','EJ301 5.51MeV Neutron recoil proton spectrum');
plot( BC501A_N5510_spectrum_exp(:,1)*EJ301_Energy_calibration(1)+EJ301_Energy_calibration(2)/1000, EJ301_N5510_spectrum_exp(:,2),'color','red','Linewidth',2);
title('EJ301 5.51MeV Neutron recoil proton spectrum','FontSize',20,'Color','black','FontWeight','bold');
xlabel('Equivalent electronic energy(MeVee)','FontSize',20,'Color','black','FontWeight','bold');
ylabel('Counts','FontSize',20, 'Color','black','FontWeight','bold');
