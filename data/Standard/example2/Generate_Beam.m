%=========================================================
%           Generate the nodes and elements in beam
%                  WangShuai 2020/4/23
%=========================================================

clear
clc

fid = fopen('Beam_data.txt','w');
for jj = 1 : 11
    fprintf(fid,'%d\t %d\t %d\t %d\t',jj, 1, 1, 1);
    fprintf(fid,'%d\t %d\t %d\t',0.0d0, (jj - 1) * 0.05, 0.0d0);
    fprintf(fid,'\r\n');
end
for ii = 2 : 41
    for jj =  1 : 11
        fprintf(fid,'%d\t %d\t %d\t %d\t',(ii - 1) * 11 + jj, 0, 0, 1);
        fprintf(fid,'%d\t %d\t %d\t',(ii - 1) * 0.05, (jj - 1) * 0.05, 0.0d0);
        fprintf(fid,'\r\n');
    end
end
 fprintf(fid,'\r\n');
for ii = 1 : 40
    for jj = 1 : 10
       fprintf(fid,'%d\t %d\t %d\t %d\t %d\t %d\t',...
           (ii-1)*10+jj, (ii-1)*11+jj, (ii-1)*11+jj+1, (ii-1)*11+jj+12, (ii-1)*11+jj+11, 1);
       fprintf(fid,'\r\n');
    end
end
 fprintf(fid,'\r\n');

fclose(fid);