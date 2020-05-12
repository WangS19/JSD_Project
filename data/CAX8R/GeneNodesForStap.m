function [Nodes,Eles] = GeneNodesForStap()
% �ڵ���Ϣ  ��ţ�1��  Լ����3��  λ�ã�3��
% ��Ԫ��Ϣ  ��Ԫ��ţ�1�� ��Ԫ�ڵ��ţ�8�� 
%% ����1/4�����Գƽṹ�Ķ�άCAX8R����
NodeNumPerLineInEle = 3;
RadEleNum = 10;           
TanEleNum = 20;
Rmin = 10;                 %�ڰ뾶
Rmax = 20;                 %��뾶
RadNum = 2*RadEleNum+1;    %����3�ڵ�����
RadNum2 = RadEleNum+1;     %����2�ڵ�����
NodesNum = (RadNum+RadNum2)*TanEleNum+RadNum;
EleNum = RadEleNum*TanEleNum;

Eles = zeros(EleNum,9);
Nodes = zeros(NodesNum,7);
L = Rmax-Rmin;
Theta = pi/2;
dL = L/2/RadEleNum;
dTheta = Theta/2/TanEleNum;

%% ��Ԫ�ڵ���Ϣ����
%�ڵ� �����Ⱦ�������������
NodeCount = 0;
for i = 1:TanEleNum
    % theta_0
    theta_0 = (i-1)*2*dTheta;
    % ����3�ڵ��һ��
    for j = 1:RadNum
        NodeCount = NodeCount+1;
        Nodes(NodeCount,1) = NodeCount;
        % Ĭ�����ɶ�3Լ����Ϊ0
        Nodes(NodeCount,4) = 1;
        Nodes(NodeCount,7) = 0;
        % ����r_p��theta_p
        theta_p = theta_0;
        r_p = Rmin +(j-1)*dL;
        % ����xp, yp
        xp = r_p*cos(theta_p);
        yp = r_p*sin(theta_p);
        Nodes(NodeCount,5) = xp;
        Nodes(NodeCount,6) = yp;
        % Լ��theta=0�����нڵ�
        if  (i == 1)
            Nodes(NodeCount,2) = 1;
            Nodes(NodeCount,3) = 1;
        end
    end
    % ����2�ڵ��һ��
    for j = 1:RadEleNum+1
        NodeCount = NodeCount+1;
        Nodes(NodeCount,1) = NodeCount;
        % Ĭ�����ɶ�3Լ����Ϊ0
        Nodes(NodeCount,4) = 1;
        Nodes(NodeCount,7) = 0;
        % ����r_p��theta_p
        theta_p = theta_0 + dTheta;
        r_p = Rmin +(j-1)*2*dL;
        % ����xp, yp
        xp = r_p*cos(theta_p);
        yp = r_p*sin(theta_p);
        Nodes(NodeCount,5) = xp;
        Nodes(NodeCount,6) = yp;
    end
    if (i == TanEleNum)
        for j = 1:RadNum
            NodeCount = NodeCount+1;
            Nodes(NodeCount,1) = NodeCount;
            % Ĭ�����ɶ�3Լ����Ϊ0
            Nodes(NodeCount,4) = 1;
            Nodes(NodeCount,7) = 0;
            % ����r_p��theta_p
            theta_p = theta_0+2*dTheta;
            r_p = Rmin +(j-1)*dL;
            % ����xp, yp
            xp = r_p*cos(theta_p);
            yp = r_p*sin(theta_p);
            Nodes(NodeCount,5) = xp;
            Nodes(NodeCount,6) = yp;
            % Լ��theta=pi/4�����нڵ��x����
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

%��Ԫ��Ϣ        
EleCount = 0;
for i = 1:TanEleNum              %����
    for j = 1:RadEleNum          %����
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