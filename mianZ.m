%%%%%%%%%%%%%GRGMF,5-fold
clc
clear
% % %%%%%%%%%%%%step1�� pair prediction
cv_flag = 1;
V_D = process_data( );
prop = 1;
main_cv(cv_flag,V_D,prop);
% 



% %%%%%%%%%%%%step2�� new virus
cv_flag = 2;
V_D = process_data( );
prop = 100;
main_cv(cv_flag,V_D,prop);

%


% %%%%%%%%%%%%step2��new drug
cv_flag = 3;
V_D = process_data( );
prop = 100;
main_cv(cv_flag,V_D,prop);
% 

