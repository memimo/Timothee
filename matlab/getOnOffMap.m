function map=getOnOffMap(idx)
global COMMON
tmp = COMMON.onOffMap{idx};

onMap = max(0,tmp);
offMap = - min(0,tmp);

tmp2(:,:,1) = onMap;
tmp2(:,:,2) = offMap;

map{1} = tmp2;
