clc
clear all
% Obtain the input data
% Water Frequency map is derived from multi-temporal RS imagery.
% The starting point (i.e., the watershed outlet or most downstream point) and ending point (i.e., the most upstream point) of a river path can be manually specified according to high-resolution Google earth imagery or the derived water frequency map.
[num,R]=geotiffread('Water_Frequency.tif'); % Reading the Water-occurrence frequency image. Here 'num' is the Water-occurrence frequency matrix and R is the spatial reference object.
[Yt,~]=importdata('Upstream_point.txt');  % Reading row and column number of most upstream points for each river within the Water-occurrence frequency image. The data in 'Upstream_point.txt' is a 2-column numeric array, i.e.,  the row and column number of each upstream point within the Water-occurrence frequency image. Each row represents an upstream point. Yt(:,1) is row number and Yt(:,2) is column number.
[Ot,~]=importdata('Outlet_point.txt');  % Reading row and column number of outlet point within the Water-occurrence frequency image. Ot(1,1) is row number and Ot(1,2) is column number.

% Set constant value
High_cost=100000000;  % Set a very high cost score value for the pixel without water.
Mim_value = 5; % Set a threshold value for water-occurrence frequency. Pixels with water-occurrence frequency less than this threshold value will be set nodata (-9999), i.e., no water occurrence in order to avoid the impact of abnormal interpretation results.
Index_a=[-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1]; % Represent 8 directions for each pixel

% Pre-process the input data
num(num<=Mim_value)=-9999;  % Find pixels with water-occurrence frequency less than threshold value and set the frequency as -9999, i.e., no water occurrence
a1=Ot(1,1);  % Obtain row number of outlet point from the input data
b1=Ot(1,2); % Obtain column number of outlet point from the input data
[r1,c1]=size(num);

% Set cost score value for the pixels without water occurrence according to the condition of their neighboring pixels
% The cost value of non-water pixels depends on the distance to the nearest water pixel.
% Cost is the cost value matrix. Cost1, Cost2, and Num1 are intermediate variables.
Cost=num;
Cost(num<=0)=0;
Cost(num>0)=1;
Cost1=num;
for mm=1:100
    ZZ{1}=[Cost(:,2:c1) zeros(r1,1)];
    ZZ{2}=[zeros(r1,1) Cost(:,1:c1-1)];
    ZZ{3}=[Cost(2:r1,:); zeros(1,c1)];
    ZZ{4}=[zeros(1,c1);Cost(1:r1-1,:)];
    ZZ{5}=[[Cost(2:r1,2:c1) zeros(r1-1,1)]; zeros(1,c1)];
    ZZ{6}=[zeros(1,c1); [Cost(1:r1-1,2:c1) zeros(r1-1,1)]];
    ZZ{7}=[zeros(1,c1);[zeros(r1-1,1) Cost(1:r1-1,1:c1-1)]];
    ZZ{8}=[[zeros(r1-1,1) Cost(2:r1,1:c1-1)]; zeros(1,c1)];
    Cost2=ZZ{1}+ZZ{2}+ZZ{3}+ZZ{4}+ZZ{5}+ZZ{6}+ZZ{7}+ZZ{8};
    zero_index=find((Cost2>0) &(Cost==0));
    Cost1(zero_index)=mm*High_cost;
    Cost(zero_index)=1;
end
Cost=Cost1;
num1=num;
% Update cost score value for the pixels without water occurrence but adjacent to water pixel
for k=1:1000
    a=0;
    for i=2:r1-1
        for j=2:c1-1
            if num1(i,j)>=0
                for mm=1:8
                    if num1(i+Index_a(mm,1),j+Index_a(mm,2))<0
                        num1(i+Index_a(mm,1),j+Index_a(mm,2))=High_cost;
                        a=1;
                    end
                end
            end
        end
    end
    if a==0
        break
    end
end
Cost(num1==High_cost)=High_cost;

% Set cost score value for the pixels with water according to their water-occurrence frequency
% Please refer to the section 2.2
Max_frequency=max(max(num));
num(num>0)=(Max_frequency-num(num>0)).^3+1;
Cost(Cost<Max_frequency+2)=num(Cost<Max_frequency+2);

% Please refer to steps 1)-8) in the section 2.3
% Find a best path to represent the river network through the cost image
% CostAcc is the minimum cumulative cost value. CostAcc_Initial, CostAcc_indexR and CostAcc_indexC are intermediate variables.
CostAcc=zeros(r1,c1);
CostAcc(a1,b1)=1;
CostAcc_Initial=zeros(r1,c1);
for mm=1:100000
    % Compare the current CostAcc image and the previous CostAcc image and identify those corresponding pixels with different values.
    [CostAcc_indexR,CostAcc_indexC]=find(CostAcc~=CostAcc_Initial);
    CostAcc_indexR(CostAcc_indexR==1)=a1;
    CostAcc_indexC(CostAcc_indexR==1)=b1;
    CostAcc_indexR(CostAcc_indexR==r1)=a1;
    CostAcc_indexC(CostAcc_indexR==r1)=b1;
    CostAcc_indexR(CostAcc_indexC==1)=a1;
    CostAcc_indexC(CostAcc_indexC==1)=b1;
    CostAcc_indexR(CostAcc_indexC==c1)=a1;
    CostAcc_indexC(CostAcc_indexC==c1)=b1;
    [rr1,cc1]=size(CostAcc_indexR);
    aaaaa=0;
    CostAcc_Initial=CostAcc;
    % For each identified pixel, examine each of its immediate neighboring pixels and calculate the cumulative cost.
    for i=1:rr1
        for j=1:8
            a=Index_a(j,1);
            b=Index_a(j,2);
            if (Cost(CostAcc_indexR(i)+a,CostAcc_indexC(i)+b)>0) &&  (( CostAcc(CostAcc_indexR(i)+a,CostAcc_indexC(i)+b)>(CostAcc(CostAcc_indexR(i),CostAcc_indexC(i))+Cost(CostAcc_indexR(i)+a,CostAcc_indexC(i)+b)*(1+abs(a)*abs(b)*0.414214))) || (CostAcc(CostAcc_indexR(i)+a,CostAcc_indexC(i)+b)==0))
                CostAcc(CostAcc_indexR(i)+a,CostAcc_indexC(i)+b)=CostAcc(CostAcc_indexR(i),CostAcc_indexC(i))+Cost(CostAcc_indexR(i)+a,CostAcc_indexC(i)+b)*(1+abs(a)*abs(b)*0.414214);
                aaaaa=2;
            end
        end
    end
    if aaaaa==0
        break
    end
end

% Please refer to the section 2.3 
% Trace back through the minimum cumulative cost imagery from goal point to starting point to extract the minimum-cost path, i.e. the river path.
CostAcc(CostAcc<0)=0;
[Num_River,~]=size(Yt);
River=zeros(r1,c1);
for i=1:Num_River
    River(Yt(i,1),Yt(i,2))=1;
end
for i=1:Num_River % Loop for each river
    a_tem=Yt(i,1);
    b_tem=Yt(i,2);
    for mm=1:18000
        a=0;
        suiji=randperm(8);
        for j=1:8
            aaa=Index_a(suiji(1,j),1);
            bbb=Index_a(suiji(1,j),2);
            if  CostAcc(a_tem+aaa,b_tem+bbb)>0 && River(a_tem+aaa,b_tem+bbb)==0 && abs(CostAcc(a_tem+aaa,b_tem+bbb)-CostAcc(a_tem,b_tem)+Cost(a_tem,b_tem)*(1+abs(aaa)*abs(bbb)*0.414214))<1
                River(a_tem+aaa,b_tem+bbb)=i;
                a_tem=a_tem+aaa;
                b_tem=b_tem+bbb;
                a=a+1;
                if a>0
                    break
                end
            end
        end
        if a==0
            break
        end
    end
end

% Output river extraction result
info = geotiffinfo('Water_Frequency.tif');
geotiffwrite('River.tif',River,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
