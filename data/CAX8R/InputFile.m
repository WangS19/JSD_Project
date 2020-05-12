clc
clear
close all
% ���ýڵ��CAX8R��Ԫ
[Nodes,Eles] = GeneNodesForStap();%��ͬ������ͬ����Ҫ����
%% Stappp CAX8R��Ԫ �����ļ����
FileName = 'CAX8R_SPHERE'
Input = fopen([FileName,'.dat'],'w');
% ������
fprintf(Input,[FileName,'\n']);
% ������͵�
NodeNum = size(Nodes,1);
fprintf(Input,'%d 1 1 1\n',NodeNum);

% �ڵ���Ϣ  ��ţ�1��  Լ����3��  λ�ã�3��
for i=1:size(Nodes,1)
    for j=1:size(Nodes,2)
        if j< 5
            fprintf(Input,'%d ',Nodes(i,j));
        else
            fprintf(Input,'%3.4f ',Nodes(i,j));
        end
        
    end
    fprintf(Input,'\n',Nodes(i,j));
end

% �غ���Ϣ  �ڵ��ţ�1��  ���ɶȱ�ţ�1��  ��ֵ��1��
%GeneLoadsForStap(Loads);%��ͬ������ͬ����Ҫ����
Loads = [size(Nodes,1) 2 -1000.0];
LoadNum = size(Loads,1);
fprintf(Input,'1 %d\n',LoadNum);
for i = 1:LoadNum
    for j = 1:size(Loads,2)
        fprintf(Input,'%d ',Loads(i,j));
    end
    fprintf(Input,'\n');
end

% ��Ԫ��Ϣ ��Ԫ��ţ�1�� ��Ԫ�ڵ��ţ�8�� ���ϱ�ţ�1��
EleType = 8;   EleNum = size(Eles,1);   MatNum =1;  GuassNum = 2;
EleMat = ones(EleNum, 1);%Ĭ��ȫ��Ϊ1�����ж��ֲ��Ͽ������������
Elastic = '2.1e11 ';    mu = '0.3 ';     rou = '7900';
fprintf(Input,'%d %d %d %d\n',EleType,EleNum,MatNum,GuassNum);
for i = 1:MatNum
    fprintf(Input,['%d ',Elastic,mu,rou,'\n'],MatNum);
end
for i = 1:EleNum
    for j = 1:size(Eles,2)
       fprintf(Input,'%d ',Eles(i,j)); 
    end
    fprintf(Input,'%d ', EleMat(i));
    fprintf(Input,'\n');
end
fprintf(Input,'stop');
%% ABAQUS CAX8R��Ԫ �����ļ����
InputAba = fopen([FileName,'.inp'],'w');
fprintf(InputAba,'*Node\n')
for i=1:size(Nodes,1)
    fprintf(InputAba,'%d, ',Nodes(i,1));
    fprintf(InputAba,'%3.4f, ',Nodes(i,5));
    fprintf(InputAba,'%3.4f\n ',Nodes(i,6));
end
fprintf(InputAba,'*Element, type=CAX8R\n')
for i=1:size(Eles,1)
    for j=1:8
        fprintf(InputAba,'%d, ',Eles(i,j));
    end
    fprintf(InputAba,'%d\n ',Eles(i,9));
end






