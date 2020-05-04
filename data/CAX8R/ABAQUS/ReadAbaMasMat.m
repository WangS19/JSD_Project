clc
clear
% ������
MassFile = importdata('CAX8R_dyn_outM_MASS2.mtx',',',0);
Mass = zeros(16,16);
for i=1:size(MassFile,1)
    row = (MassFile(i,1)-1)*2+MassFile(i,2);
    col = (MassFile(i,3)-1)*2+MassFile(i,4);
    Mass(row,col) = MassFile(i,5);
    Mass(col,row) = Mass(row,col);
end
%����������
MassSum = sum(Mass,2);
MassTol = sum(MassSum)/2;
%�ж�����
eigMass = eig(Mass)
%�ж������غ�
MassTol-2*pi*50*7900*10000
%�Խ�������
MassDig = 0;
for i = 1:size(Mass,1)
    MassDig = MassDig+Mass(i,i);
end
MassDig/2
