%*********************
% COMPUTE ON-OFF MAP *
%*********************

for i=1:PARAM.c.RFSize(1)
    internalPos(1) = 1 + (i-floor((PARAM.c.RFSize(1)+1)/2))*PARAM.s.subSampling*PARAM.onOff.subSampling;
    for j=1:PARAM.c.RFSize(2)
        timedLogLn(['ON-OFF MAP (' int2str(i) ',' int2str(j) ')']);
        internalPos(2) = 1 + (j-floor((PARAM.c.RFSize(2)+1)/2))*PARAM.s.subSampling*PARAM.onOff.subSampling;
        onOffMap = computeOnOffMap(PARAM.DoGFilter,PATH.image,PARAM.cropPos,internalPos,PARAM.s.RFSize,PARAM.onOff.subSampling,PARAM.zoom);
        save([PATH.map 'onOff.' int2str(i) '.' int2str(j) '.mat'],'onOffMap')
        clear onOffMap
    end
end

