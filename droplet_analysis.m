tic
clear all;
disp ('��͊J�n')

Number=0;
nimus_NUM=0;
pix=1; %640�Ȃ�100/85�@�}�C�N��/�s�N�Z���@532��1 ���� 640��1�ōs����

bot_num=1;

DIRECTORY=pwd;
address=DIRECTORY;
tifFiles=dir('*.tif');
picture_num=numel(tifFiles);
allowable_1=30;
allowable_2=30;
allowable_3=30;
Number_M=0;
a=1;
chinp=1;
% PON_THR=zeros(255,3);
%for a=1%:picture_num
% mkdir('ryuukei');
%     mkdir('�������]');
%     mkdir('���');
%     mkdir('��͂Q');

data = imread((tifFiles(a).name));
disp(tifFiles(a).name);
siguma=300;
sz = size(data);

% Diameter=zeros(5000000,8);
% BOTSU=zeros(5000000,8);

atta=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%������臒l�X�^�[�g�ݒ�
for PON= 1:170
    data1 = imbinarize(data,  PON/255); %��l��
    data1 = imcomplement(data1); %�␔�@�������]
    data1 = imfill(data1,'holes');%�������ӏ���⊮
    
    %     imwrite(data1,strcat(char(address),'/��l��/',tifFiles(a).name),'TIFF');
    data3 = imcomplement(data);
    %     imwrite(data3,strcat(char(address),'/�������]/',tifFiles(a).name),'TIFF'); %�������݂��Ă�
    
    CC(PON) = bwconncomp(data1, 8); %�Ȃɂ��Ă�̂��� �񎟌��̉摜�ł��I���Č����Ă�
    
    Characteristic  = regionprops(CC(PON),'Centroid','Circularity','EquivDiameter');
    %'Area'�̈���̎��ۂ̃s�N�Z����
    %'Centroid'�̈�̏d�S
    %'Circularity'�I�u�W�F�N�g�̊ۂ����w�肷��^�~�x(4*Area*pi)/(Perimeter^2)
    %'Eccentricity'�̈�Ɠ��� 2 �����[�����g�����ȉ~�̗��S��
    %'EquivDiameter'sqrt(4*Area/pi) �����~���a
    %'Perimeter'�̈�̋��E�̎��͒�
    
    %Characteristic����f�[�^�̒��o
    Centroid=cat(1,Characteristic.Centroid)*pix;
    Circularity=cat(1,Characteristic.Circularity);
    EquivDiameter=cat(1,Characteristic.EquivDiameter)*pix;
    %Characteristic����f�[�^�̒��o�I
    
    nC=CC(PON).NumObjects;
    for aaa = 1:CC(PON).NumObjects
        if( Circularity(aaa,1)<0.5  ||EquivDiameter(aaa,1)>500  || Centroid(aaa,2)<1000)
            nC=nC-1;
            data1(CC(PON).PixelIdxList{aaa}) = 0;
        end
    end
    
    CC(PON).NumObjects
    %nC
    
    CC(PON) = bwconncomp(data1, 8);
    
    Characteristic  = regionprops(CC(PON),'Centroid','Circularity','EquivDiameter');
    %'Area'�̈���̎��ۂ̃s�N�Z����
    %'Centroid'�̈�̏d�S
    %'Circularity'�I�u�W�F�N�g�̊ۂ����w�肷��^�~�x(4*Area*pi)/(Perimeter^2)
    %'Eccentricity'�̈�Ɠ��� 2 �����[�����g�����ȉ~�̗��S��
    %'EquivDiameter'sqrt(4*Area/pi)
    %'Perimeter'�̈�̋��E�̎��͒�
    
    %Characteristic����f�[�^�̒��o
    Centroid=cat(1,Characteristic.Centroid);
    Circularity=cat(1,Characteristic.Circularity);
    EquivDiameter=cat(1,Characteristic.EquivDiameter);
    
    for a1 = 100:10100
        for a2 = 500:600
            data1(a1,a2) = 255;
        end
    end
    
    %     imwrite(data1,strcat(char(address),'/���/',tifFiles(a).name),'TIFF');
    %%%�Ƃ肠���������܂�
    nn = 0;
    
    mukimuki=0;
    for aaa = 1:CC(PON).NumObjects %��
        nn = nn + 1;
        n_1 = 0;
        base1 = 0;
        if(Centroid(aaa,1)>500 && Centroid(aaa,2)>500 && Centroid(aaa,1)<sz(2)-500 && Centroid(aaa,2)<sz(1)-500) % if ��
            for ii = round(Centroid(aaa,1) - 3*EquivDiameter(aaa,1)/4):round(Centroid(aaa,1) + 3*EquivDiameter(aaa,1)/4)
                for jj = round(Centroid(aaa,2) - 3*EquivDiameter(aaa,1)/4):round(Centroid(aaa,2) + 3*EquivDiameter(aaa,1)/4)
                    if (ii - Centroid(aaa,1))^2 + (jj - Centroid(aaa,2))^2 >= ((3*EquivDiameter(aaa,1)/4)-1)^2 && (ii - Centroid(aaa,1))^2 + (jj - Centroid(aaa,2))^2 <= ((3*EquivDiameter(aaa,1)/4)+1)^2
                        n_1 = n_1 +1;
                        base1 = base1 + double(data3(jj,ii));
                    end
                end
            end
            base1 = double(base1 / n_1);
            
            n_2 = 0;
            base2 = 0;
            for iii = round(Centroid(aaa,1) - EquivDiameter(aaa,1)/4):round(Centroid(aaa,1) + EquivDiameter(aaa,1)/4)
                for jjj = round(Centroid(aaa,2) - EquivDiameter(aaa,1)/4):round(Centroid(aaa,2) + EquivDiameter(aaa,1)/4)
                    if (iii - Centroid(aaa,1))^2 + (jjj - Centroid(aaa,2))^2 >= ((EquivDiameter(aaa,1)/4)-1)^2  && (iii - Centroid(aaa,1))^2 + (jjj - Centroid(aaa,2))^2 <= ((EquivDiameter(aaa,1)/4)+1)^2
                        n_2 = n_2 +1;
                        base2 = base2 + double(data3(jjj,iii));
                    end
                end
            end
            base2 = double(base2 / n_2);
           
            
            n_4 = 0;
            II = 0;
            for i = round(Centroid(aaa,1) - 3*EquivDiameter(aaa,1)/4):round(Centroid(aaa,1) +3*EquivDiameter(aaa,1)/4)
                for j = round(Centroid(aaa,2) - 3*EquivDiameter(aaa,1)/4):round(Centroid(aaa,2) + 3*EquivDiameter(aaa,1)/4)
                    if (i - Centroid(aaa,1))^2 + (j - Centroid(aaa,2))^2 <= (3*EquivDiameter(aaa,1)/4)^2 && data3(j,i) < (base1 + base2)/2
                        n_4 = n_4 +1;
                        I = (double(data3(j,i)) - base1)/(base2 - base1);
                        %                         II = II + (I - I_ave)^2;
                        II = II + (I - 0)^2;
                    end
                end
            end
            I_var1= II / (n_4 -1);
            
            n_6 = 0;
            III = 0;
            for i = round(Centroid(aaa,1) - 3*EquivDiameter(aaa,1)/4):round(Centroid(aaa,1) + 3*EquivDiameter(aaa,1)/4)
                for j = round(Centroid(aaa,2) - 3*EquivDiameter(aaa,1)/4):round(Centroid(aaa,2) + 3*EquivDiameter(aaa,1)/4)
                    if (i - Centroid(aaa,1))^2 + (j - Centroid(aaa,2))^2 <= (3*EquivDiameter(aaa,1)/4)^2 && data3(j,i) > (base1 + base2)/2
                        n_6 = n_6 +1;
                        I = (double(data3(j,i)) - base1)/(base2 - base1);
                        %                         II = II + (I - I_ave)^2;
                        III = III + (I - 1)^2;
                    end
                end
            end
            I_var2 = III / (n_6 -1);
            
            I_var = I_var1 + I_var2 ;
            
            fd=-1*(I_var-0.505*EquivDiameter(aaa,1)^(-0.305))/(0.105*EquivDiameter(aaa,1)^(-0.093));
            
            Area=0;
            for iiii = round(Centroid(aaa,1) - EquivDiameter(aaa,1)):round(Centroid(aaa,1) + EquivDiameter(aaa,1))
                for jjjj = round(Centroid(aaa,2) - EquivDiameter(aaa,1)):round(Centroid(aaa,2) + EquivDiameter(aaa,1))
                    if (iiii - Centroid(aaa,1))^2 + (jjjj - Centroid(aaa,2))^2 <= (EquivDiameter(aaa,1))^2 && (255- PON) < data3(jjjj,iiii)
                        Area=Area+1;
                    end
                    EquivalentDiameter2=sqrt(4*Area/pi);
                end
            end
            
            if(fd>0.8 && EquivalentDiameter2/EquivDiameter(aaa,1)<1.2 && EquivalentDiameter2/EquivDiameter(aaa,1)>0.8 && base2>base1+35)
                
                if atta==0
                    Number=Number+1;
                    Diameter(Number,1)=Number;                    %�ԍ�
                    Diameter(Number,2)=EquivDiameter(aaa,1);      %���a
                    Diameter(Number,3)=Centroid(aaa,1);           %x�t��� �t���Ⴄ��
                    Diameter(Number,4)=Centroid(aaa,2);           %y
                    if fd<=1
                        OPP=log(fd);
                        XXX=-2*siguma*siguma*OPP;
                        z=sqrt(XXX);
                        Diameter(Number,5)=z;
                    else
                        Diameter(Number,5)=0;
                    end
                    Diameter(Number,6)=PON;                       %臒l
                    
%                     z=sqrt(-2*siguma*siguma*log(fd));
%                     Diameter(Number,6)=z;                         %z��Βl
                    Diameter(Number,7)=fd;                        %fd
                end
                if atta>0
                    AKAISHI=0;
                    KOHEI=0;
                    MUNYA=0;
%                     Number=PON_THR(atta,2); %������������
                    for Redstone=1 : Number %mom
                        %                     Redstone
                        if (Centroid(aaa,1)-Diameter(Redstone,3))^2+(Centroid(aaa,2)-Diameter(Redstone,4))^2<(EquivDiameter(aaa,1)/2)^2 %) && (abs(Centroid(aaa,2)-Diameter(Redstone,3))<allowable_3)) %%%%%%%%%%if�̒��ɗ��a���f���������Ǐ�����abs(EquivDiameter(aaa,1)-Diameter(Redstone,1))<allowable_1 &&
                            AKAISHI=AKAISHI+1; %�L���̒��ɂ�������A�J�C�V�P
                            KOHEI=1;
                            MUNYA=1;
                            memo=Redstone; %memo��Diameter�ɂ���o�g���i���o�[
%                             break
                            
                            if AKAISHI==1 && MUNYA==1 %AKAISHI���P�@�L���ɂ������D
                                %                         BOTSU(bot_num,1)=Diameter(memo,1);
                                %                         BOTSU(bot_num,2)=Diameter(memo,2);
                                %                         BOTSU(bot_num,3)=Diameter(memo,3); %BOTSUhairetsu_redstone.m
                                %                         BOTSU(bot_num,4)=Diameter(memo,4);
                                %                         BOTSU(bot_num,5)=memo;
                                %                         BOTSU(bot_num,6)=Diameter(memo,6);
                                %                         BOTSU(bot_num,7)=Diameter(memo,7);
                                %                         BOTSU(bot_num,8)=bot_num;
                                %���ꂩ��
                                Diameter(memo,2)=EquivDiameter(aaa,1);
                                Diameter(memo,3)=Centroid(aaa,1);
                                Diameter(memo,4)=Centroid(aaa,2);
                                if fd<=1
                                    OPP=log(fd);
                                    XXX=-2*siguma*siguma*OPP;
                                    z=sqrt(XXX);
                                    Diameter(memo,5)=z;
                                else
                                    Diameter(memo,5)=0;
                                end
                                Diameter(memo,6)=PON;
                                %                         z=sqrt(-2*siguma*siguma*log(fd));
                                %                         Diameter(Number,6)=z;
                                Diameter(memo,7)=fd;
                                %                         Diameter(memo,8)=bot_num;
                                bot_num=bot_num+1;
                                MUNYA=0;
                                
                            elseif AKAISHI>1 && MUNYA==1 %AKAISHI���P�@�L���ɂ������D
                                %-100�����
                                Diameter(memo,2)=0;%EquivDiameter(aaa,1);
                                Diameter(memo,3)=0;%Centroid(aaa,1);
                                Diameter(memo,4)=0;%Centroid(aaa,2);
%                                 if fd<=1
%                                     OPP=log(fd);
%                                     XXX=-2*siguma*siguma*OPP;
%                                     z=sqrt(XXX);
                                Diameter(memo,5)=0;%z;
%                                 end
                                Diameter(memo,6)=PON;
                                %                         z=sqrt(-2*siguma*siguma*log(fd));
                                %                         Diameter(Number,6)=z;
                                Diameter(memo,7)=0;%fd;
                                %                         Diameter(memo,8)=bot_num;
                                nimus_NUM=nimus_NUM+1;
                                MUNYA=0;
                            end
                        end
                    end %mom
                    
                    if KOHEI==0
                        Number=Number+1;
                        Diameter(Number,1)=Number;
                        Diameter(Number,2)=EquivDiameter(aaa,1);
                        Diameter(Number,3)=Centroid(aaa,1);
                        Diameter(Number,4)=Centroid(aaa,2);
                        if fd<=1
                            OPP=log(fd);
                            XXX=-2*siguma*siguma*OPP;
                            z=sqrt(XXX);
                            Diameter(Number,5)=z;
                        else
                            Diameter(Number,5)=0;
                        end
                        
                        Diameter(Number,6)=PON;
                        
                        %                         z=sqrt(-2*siguma*siguma*log(fd));
                        %                         Diameter(Number,6)=z;
                        Diameter(Number,7)=fd;
                        mukimuki=mukimuki+1;
                    end
                    
                    
                    
                end
            end
        end %if ��
    end %��
    PON
    if Number>0 && atta==0
        atta=atta+1;
        %         PON_THR(atta,1)=PON;
        %         PON_THR(atta,2)=Number;
        %         PON_THR(atta,3)=mukimuki;
%         PON
        Number
        %         disp('KOKO');
    end
    if mukimuki>0
        atta=atta+1;
        %         PON_THR(atta,1)=PON;
        %         PON_THR(atta,2)=Number;
        %         PON_THR(atta,3)=mukimuki;
%         PON
        Number
        %         disp('KOKO111111');
    end
end
nimus_NUM


writematrix(Diameter,'001_adopt_Diameter532.txt');

Histgram=zeros(500,2);
for i=1:Number
    for j=1:500
        if(j==round(Diameter(i,2))&&Diameter(i,2)>0)
            Histgram(j,2)=Histgram(j,2)+1;
        end
    end
end
for j=1:500
    Histgram(j,1)=j;
end
d1 = 0;
for i1 = 1:500
    d1 = d1 + Histgram(i1,2)*i1;
end

Number=Number-nimus_NUM;
Number

d10 = d1/Number;

d3 = 0;
d2 = 0;
for i2 = 1:500
    d3 = d3 + Histgram(i2,2)*i2*i2*i2;
    d2 = d2 + Histgram(i2,2)*i2*i2;
end
d32 = d3/d2;
%end

HOZON=zeros(4,1);
HOZON(1,1)=d10;
HOZON(2,1)=d32;
HOZON(3,1)=Number;
HOZON(4,1)=nimus_NUM;
writematrix(HOZON,'005_d10_d32_Number532_minus.txt');
d10
d32

Q_spray=0;
for i3 = 1:500
    Q_spray = Q_spray+Histgram(i3,2)*i3*10^(-6)*i3*10^(-6)*i3*10^(-6)*(1/6)*695000*pi;
end

bot_num=bot_num-1;


% writematrix(PON_THR,'003_threshold.txt');
writematrix(Histgram,'004_HISTGRAM532.txt');

data4 = zeros(round(sz(1)/10)+1500,round(sz(2)/10));

for i = 1:round(sz(1)/10)
    for j = 1:round(sz(2)/10)
        for ii=1:Number
            if Diameter(ii,2)>0
                if (i - round(Diameter(ii,4)/10))^2 + (j - round(Diameter(ii,3)/10))^2 < 500
                    data4(i,j) =Diameter(ii,2);
                end
            end
        end
    end
    if(rem(i,500)==0)
        i
    end
disp('�I')
toc
clear all;