clc
clear
close all
% 设置节点和CAX8R单元
[Nodes,Eles] = GeneNodesForStap();%不同算例不同，需要更改
%% Stappp CAX8R单元 输入文件输出
FileName = 'CAX8R_SPHERE'
Input = fopen([FileName,'.dat'],'w');
% 标题行
fprintf(Input,[FileName,'\n']);
% 求解类型等
NodeNum = size(Nodes,1);
fprintf(Input,'%d 1 1 1\n',NodeNum);

% 节点信息  编号（1）  约束（3）  位置（3）
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

% 载荷信息  节点编号（1）  自由度编号（1）  数值（1）
%GeneLoadsForStap(Loads);%不同算例不同，需要更改
Loads = [size(Nodes,1) 2 -1000.0];
LoadNum = size(Loads,1);
fprintf(Input,'1 %d\n',LoadNum);
for i = 1:LoadNum
    for j = 1:size(Loads,2)
        fprintf(Input,'%d ',Loads(i,j));
    end
    fprintf(Input,'\n');
end

% 单元信息 单元编号（1） 单元节点标号（8） 材料编号（1）
EleType = 8;   EleNum = size(Eles,1);   MatNum =1;  GuassNum = 2;
EleMat = ones(EleNum, 1);%默认全部为1，含有多种材料可以针对性设置
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
%% ABAQUS CAX8R单元 输入文件输出
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






