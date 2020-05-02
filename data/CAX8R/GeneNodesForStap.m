function [Nodes,Eles] = GeneNodesForStap()
% 节点信息  编号（1）  约束（3）  位置（3）
% 单元信息  单元编号（1） 单元节点标号（8） 
%% 建立1/4球壳轴对称结构的二维CAX8R网格
NodeNumPerLineInEle = 3;
RadEleNum = 10;           
TanEleNum = 20;
Rmin = 10;                 %内半径
Rmax = 20;                 %外半径
RadNum = 2*RadEleNum+1;    %径向3节点数量
RadNum2 = RadEleNum+1;     %径向2节点数量
NodesNum = (RadNum+RadNum2)*TanEleNum+RadNum;
EleNum = RadEleNum*TanEleNum;

Eles = zeros(EleNum,9);
Nodes = zeros(NodesNum,7);
L = Rmax-Rmin;
Theta = pi/2;
dL = L/2/RadEleNum;
dTheta = Theta/2/TanEleNum;

%% 单元节点信息计算
%节点 按照先径向再切向排列
NodeCount = 0;
for i = 1:TanEleNum
    % theta_0
    theta_0 = (i-1)*2*dTheta;
    % 含有3节点的一行
    for j = 1:RadNum
        NodeCount = NodeCount+1;
        Nodes(NodeCount,1) = NodeCount;
        % 默认自由度3约束且为0
        Nodes(NodeCount,4) = 1;
        Nodes(NodeCount,7) = 0;
        % 计算r_p，theta_p
        theta_p = theta_0;
        r_p = Rmin +(j-1)*dL;
        % 计算xp, yp
        xp = r_p*cos(theta_p);
        yp = r_p*sin(theta_p);
        Nodes(NodeCount,5) = xp;
        Nodes(NodeCount,6) = yp;
        % 约束theta=0的所有节点
        if  (i == 1)
            Nodes(NodeCount,2) = 1;
            Nodes(NodeCount,3) = 1;
        end
    end
    % 含有2节点的一行
    for j = 1:RadEleNum+1
        NodeCount = NodeCount+1;
        Nodes(NodeCount,1) = NodeCount;
        % 默认自由度3约束且为0
        Nodes(NodeCount,4) = 1;
        Nodes(NodeCount,7) = 0;
        % 计算r_p，theta_p
        theta_p = theta_0 + dTheta;
        r_p = Rmin +(j-1)*2*dL;
        % 计算xp, yp
        xp = r_p*cos(theta_p);
        yp = r_p*sin(theta_p);
        Nodes(NodeCount,5) = xp;
        Nodes(NodeCount,6) = yp;
    end
    if (i == TanEleNum)
        for j = 1:RadNum
            NodeCount = NodeCount+1;
            Nodes(NodeCount,1) = NodeCount;
            % 默认自由度3约束且为0
            Nodes(NodeCount,4) = 1;
            Nodes(NodeCount,7) = 0;
            % 计算r_p，theta_p
            theta_p = theta_0+2*dTheta;
            r_p = Rmin +(j-1)*dL;
            % 计算xp, yp
            xp = r_p*cos(theta_p);
            yp = r_p*sin(theta_p);
            Nodes(NodeCount,5) = xp;
            Nodes(NodeCount,6) = yp;
            % 约束theta=pi/4的所有节点的x方向
            Nodes(NodeCount,2) = 1;
        end
    end
end
for i = 1:size(Nodes,1)
     plot(Nodes(i,5),Nodes(i,6),'ko'); hold on;
     if (Nodes(i,2) == 1 && Nodes(i,3) == 1 )
         plot(Nodes(i,5),Nodes(i,6),'ro'); hold on;
     end
end

%单元信息        
EleCount = 0;
for i = 1:TanEleNum              %切向
    for j = 1:RadEleNum          %径向
        EleCount = EleCount+1;
        Eles(EleCount,1) =  EleCount;
        %NodeID
        NodeID0 = (i-1)*(RadNum+RadNum2)+1+(j-1)*2;
        NodeID = zeros(8,1);
        NodeID(1) = NodeID0;
        NodeID(8) = NodeID0+(RadEleNum-j+1)*2+j;
        NodeID(4) = NodeID0+RadNum+RadNum2;
        NodeID(5) = NodeID(1)+1;
        NodeID(2) = NodeID(1)+2;
        NodeID(6) = NodeID(8)+1;
        NodeID(7) = NodeID(4)+1;
        NodeID(3) = NodeID(4)+2;
        %Ele
        for k = 1:size(NodeID,1)
            Eles(EleCount,k+1) = Nodes(NodeID(k),1);
        end
    end
end
for i=1:8
    plot( Nodes(Eles(160,i+1),5), Nodes(Eles(160,i+1),6),'go');hold on;
end