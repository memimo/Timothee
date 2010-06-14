function [weight, evol, inputThreshold, nFiringInput, nFiringOutput, nAboveThr]=learnInvariance(getMap,weight,PARAM)
% Timothee Masquelier timothee.masquelier@alum.mit.edu August 2006
% see Einhauser 2002, Spratling 2005
rand('twister',5489);

CORREL_PREV_INPUT_CURR_OUTPUT = false;  % when true, reinforce when previous input and current output correlate (eg Einhauser 2002)
% when false, reinforce when previous output and current input correlate(eg Spratling 2005)
% false much more robust (Tim 06/2007)
% global PARAM
% global COMMON

% PARAM.nu = 100;
% PARAM.thresoldAM = .05;
% PARAM.thrDecay = 2^-11; % optimal around 2^-11
geomFactor = ( PARAM.learningRate_f / PARAM.learningRate_i )^(1/floor(PARAM.nIter/1000));
PARAM.a_pos = PARAM.learningRate_i;
PARAM.a_neg = -PARAM.a_pos/(PARAM.nFeat*prod(PARAM.RFSize))*PARAM.correctRatio;

if isempty(weight)
    for n=1:PARAM.n
        % format: i x j x feat x node
        weight(:,:,:,n) = getInitialWeight(PARAM.RFSize,PARAM.nFeat);
    end
    % init middle layer
    %     trAM = 1/PARAM.nu*ones([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat
    if PARAM.inputNu < Inf
        trInput = 1/PARAM.inputNu*ones([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat
    else
        trInput = ones([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat;
    end
    if PARAM.outputNu < Inf
        trOutput = 1/PARAM.outputNu*ones(1,PARAM.n);
    else
        trOutput = ones(1,PARAM.n);
    end
    %     thresholdAM = .00*ones([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat
    %     inputThreshold = 0.055*ones([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat
    inputThreshold = zeros([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat
    outputThreshold = ones(1,PARAM.n);

    previousWinnerPair = -1*ones(1,4); % format: inputWinnerI winnerInputJ winnerInputF previousOutputWinner

    %     previousMaxAM = -1;
    %     previousWinnerMI = -1; % has no impact
    %     previousWinnerMJ = -1; % has no impact
    %     previousWinnerMF = -1; % has no impact

    % init top layer
    %     trAT = 1/PARAM.nu*ones(1,PARAM.n);
    %     previousMaxAT = -1;
    %     previousWinnerT = -1;

    %reporting
    nFiringOutput = zeros(1,PARAM.n);
    nFiringInput = zeros([PARAM.RFSize PARAM.nFeat]);
    nAboveThr = zeros(1,PARAM.nIter);

    % init evol
    evol = zeros(1,PARAM.n,floor(PARAM.nIter/100));
end
disp('Learning invariances...')
tic
offset = floor(rand*PARAM.nIter);
% disp(int2str(offset))

for i=1:PARAM.nIter

    effective = false;

    % CURRENT INPUT
    grossInput = feval(getMap,i+offset); % format: i x j x feat

    % update traces
    %     if ~PARAM.trInputEffective || effective % only compute average with effective stim
    %         trInput = input / PARAM.inputNu + (1-1/PARAM.inputNu) * trInput;
    trInput = grossInput / PARAM.inputNu + (1-1/PARAM.inputNu) * trInput;
    %     end

    % normalization by the trace
    input = grossInput ./ trInput; % format: i x j x feat

    % n above thr
    nAboveThr(i) = sum(sum(sum(input>inputThreshold)));
    %     fprintf(1,'%d ', nAboveThr(i));

    % get max input
    [ inputMax idx ] =  multiDimensionalMax(input);
    inputWinnerI = idx(1);
    inputWinnerJ = idx(2);
    inputWinnerF = idx(3);
    grossInputMax = grossInput(idx(1),idx(2),idx(3));

    % CURRENT OUTPUT
    
    for n=1:PARAM.n
        if PARAM.useSoftMax
            % softmax
            w = weight(:,:,:,n);
            grossOutput(n) = softmax(w(:),input(:),PARAM.p,PARAM.q,PARAM.r);
        else
            grossOutput(n) = max(max(max( weight(:,:,:,n).*input )));
        end
    end

    % update traces
    %     if ~PARAM.trOutputEffective || effective % only compute average with effective stim
    %         trOutput = output / PARAM.outputNu + (1-1/PARAM.outputNu) * trOutput;
    if i==1
        trOutput = grossOutput;
    else
        trOutput = grossOutput / PARAM.outputNu + (1-1/PARAM.outputNu) * trOutput;
    end
    %     end

    % NOT USED ANYMORE update output thr
%     outputThreshold(previousOutputWinner) = outputThreshold(previousOutputWinner)*(previousOutputMax/PARAM.outputNu + (1-1/PARAM.outputNu) ) ;

    % normalization by the trace
    for n=1:PARAM.n
        output(n) = grossOutput(n) / trOutput(n);
    end

    % find winner in top layer
%     [outputMax outputWinner] = max(output);
            [outputMax outputWinner] = max(trOutput);

    if CORREL_PREV_INPUT_CURR_OUTPUT % correlation between previous input and current output (like Einhauser 2002)

        % learning
        %         if i>1 && previousGrossInputMax > inputThreshold(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF)
        if i>1 && previousInputMax > inputThreshold(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF)
            %             effective = true;
            % update thr
            inputThreshold(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF) = ...
                PARAM.thrCoef*previousInputMax;
            %                 PARAM.thrCoef*previousGrossInputMax;
            % update weights:
            % store winner-winner weight
            winWeight = weight(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF,outputWinner);
            
            % depress everyone
            weight(:,:,:,outputWinner) = weight(:,:,:,outputWinner) + PARAM.a_neg*weight(:,:,:,outputWinner).*(1-weight(:,:,:,outputWinner));
            % potentiate synapse between winners
            weight(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF,outputWinner) = winWeight + PARAM.a_pos*winWeight*(1-winWeight);
            
%             % Einhauser 2002
%             % depress everyone
%             weight(:,:,:,outputWinner) = weight(:,:,:,outputWinner) - PARAM.a_pos/256*weight(:,:,:,outputWinner);
%             % potentiate synapse between winners
%             weight(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF,outputWinner) = winWeight + PARAM.a_pos/256*(1-winWeight);

            
            % reporting
            nFiringInput(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF) = ...
                nFiringInput(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF)+ 1;
            nFiringOutput(outputWinner) = nFiringOutput(outputWinner) + 1;
        end % learning

    else % correlation between previous output and current input

        % learning
        if i>1 && inputMax > inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF)
            %         if i>1 && grossInputMax > inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF)
            %             effective = true;
            % update thr
            inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF) = PARAM.thrCoef*inputMax;
            %             inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF) = PARAM.thrCoef*grossInputMax;
            %             inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF) = inputMax * inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF);
            %             inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF) = inputMax;
            % update weights:
            %             if sum( [inputWinnerI inputWinnerJ inputWinnerF previousOutputWinner] ~= previousWinnerPair ) ~= 0 % different pair of winners
            % store winner-winner weight
            winWeight = weight(inputWinnerI,inputWinnerJ,inputWinnerF,previousOutputWinner);

            % depress everyone
            weight(:,:,:,previousOutputWinner) = weight(:,:,:,previousOutputWinner) + PARAM.a_neg*weight(:,:,:,previousOutputWinner).*(1-weight(:,:,:,previousOutputWinner));
            % potentiate synapse between winners
            weight(inputWinnerI,inputWinnerJ,inputWinnerF,previousOutputWinner) = winWeight + PARAM.a_pos*winWeight*(1-winWeight);


%             % Foldiak 1991
%             weight(:,:,:,previousOutputWinner) = weight(:,:,:,previousOutputWinner) + PARAM.a_pos * 2^-8 * ...
%                         (input-weight(:,:,:,previousOutputWinner));
            
            % reporting
            nFiringInput(inputWinnerI,inputWinnerJ,inputWinnerF) = ...
                nFiringInput(inputWinnerI,inputWinnerJ,inputWinnerF)+ 1;
            nFiringOutput(previousOutputWinner) = nFiringOutput(previousOutputWinner) + 1;

            % store winner for future comparison
            previousWinnerPair = [inputWinnerI inputWinnerJ inputWinnerF previousOutputWinner];
            %             end % different winner
        end % learning


    end % switch between the 2 possible correlations
    
    weight = max(weight,0); % just in case
    weight = min(weight,1); % just in case

    %     for ti=1:size(grossInput(1))% i x j x feat
    %         for tj=1:size(grossInput(2))
    %             for tf=1:size(grossInput(3))
    %                 if input(ti,tj,tf)>inputThreshold(ti,tj,tf)
    % %                 if grossInput(ti,tj,tf)>inputThreshold(ti,tj,tf)
    %                     inputThreshold(ti,tj,tf) = PARAM.thrCoef*input(ti,tj,tf);
    % %                     inputThreshold(ti,tj,tf) = PARAM.thrCoef*grossInput(ti,tj,tf);
    %                 end
    %             end
    %         end
    %     end

    if mod(i,100)==0
        %         evol(:,:,i/100) = reshape(sum(sum(weight,1),2),[PARAM.nFeat PARAM.n]); % format: feat x node x iter
        evol(1,:,i/100) = sum(sum(sum(weight,1),2),3); % format: 1 x node x iter
    end

    %     % thr decay
    inputThreshold = inputThreshold*(1-PARAM.thrDecay);
    
    % previous := current
    previousInputMax = inputMax;
    previousInputWinnerI = inputWinnerI;
    previousInputWinnerJ = inputWinnerJ;
    previousInputWinnerF = inputWinnerF;
    previousOutputWinner = outputWinner;
    previousOutputMax = outputMax;
    previousGrossInputMax = grossInputMax;

    if mod(i,1000)==0
        PARAM.a_pos = PARAM.a_pos * geomFactor;
        PARAM.a_neg = PARAM.a_neg * geomFactor;
    end
    if mod(i,round(PARAM.nIter/100))==0
        fprintf(1,'.');
    end
end

fprintf(1,'\n');
toc
%debug
disp(['tr ' num2str(trOutput)]);
disp(['thr ' num2str(outputThreshold)]);
disp(['out ' num2str(output)]);

function weight = getInitialWeight(patchSize,nFeat)
% weight = 1 + .01*rand([patchSize nFeat]);
% weight = max(0,weight);
% weight = weight / sum(weight(:));
weight = .75*( 1 + .0 * (rand([patchSize nFeat])-.5) );
% weight = .5 * ones([patchSize nFeat]);

% function input = tmp(patchSize,nFeat)
% input = rand(patchSize,patchSize,nFeat);
