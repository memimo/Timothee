function map=getS1Map(idx)
global COMMON
global PARAM
global PATH

globalIdx = mod(idx-1,COMMON.nMap)+1; % in case nIter > nMap
block = ceil( globalIdx / PARAM.c.block ); % block number

% load new map if necessary
if ~isfield(COMMON,'currentBlock') || block ~= COMMON.currentBlock
%     disp(['Loading new map. globalIdx=' int2str(globalIdx) '. block=' int2str(block)])
%     if isfield(COMMON,'currentBlock')
%         disp(['(former block=' int2str(COMMON.currentBlock) ')'])
%     end
    load([PATH.map 's1.' sprintf('%02d',block) '.' PARAM.s.type '.mat']);
    COMMON.s1Map = s1Map;
    COMMON.currentBlock = block;
%     COMMON.shuffle = randperm(PARAM.c.block);
end

idxInBlock = globalIdx - (block-1)*PARAM.c.block;

map = COMMON.s1Map(:,:,:,idxInBlock);
% map = COMMON.s1Map(:,:,:,COMMON.shuffle(idxInBlock));
