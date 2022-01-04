function [SRvec,SMBvec,Hvec,Tsvec,Vvec,Thetavec] = readData_AlongShearMargin(icestream)

if strcmp(icestream,'bindschadler')
    SRvec = load('bind_southernmargin_strainrates.csv');
    SMBvec = load('bind_southernmargin_smb.csv');
    Hvec = load('bind_southernmargin_thickness.csv');
    Tsvec = load('bind_southernmargin_surftemp.csv');
    Vvec = load('bind_southernmargin_speed.csv');
elseif strcmp(icestream,'byrd')
    SRvec = load('byrd_northernmargin_strainrates.csv');
    SMBvec = load('byrd_northernmargin_smb.csv');
    Hvec = load('byrd_northernmargin_thickness.csv');
    Tsvec = load('byrd_northernmargin_surftemp.csv');
    Vvec = load('byrd_northernmargin_speed.csv');
elseif strcmp(icestream,'pig')
    SRvec = load('pig_southernmargin_strainrates.csv');
    SMBvec = load('pig_southernmargin_smb.csv');
    Hvec = load('pig_southernmargin_thickness.csv');
    Tsvec = load('pig_southernmargin_surftemp.csv');
    Vvec = load('pig_southernmargin_speed.csv');
else
    print('No ice stream with that code')
end

SRvec = abs(SRvec)./3.15e7;
SRvec(SRvec>1e-5) = 1e-20;
Hvec(Hvec<10) = 10;
Tsvec(Tsvec>273) = 200; 
Vvec(Vvec<0) = 0;
 
end

