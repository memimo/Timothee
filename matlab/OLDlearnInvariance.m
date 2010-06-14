function [weight, evol, inputThreshold, nFiringInput, nFiringOutput]=learnInvariance(getMap,weight,PARAM)
% Timothee Masquelier timothee.masquelier@alum.mit.edu August 2006
% see Einhauser 2002, Spratling 2005

CORREL_PREV_INPUT_CURR_OUTPUT = false;  % when true, reinforce when previous input and current output correlate (eg Einhauser 2002)
                                        % when false, reinforce when previous output and current input correlate(eg Spratling 2005)
% global PARAM
% global COMMON

% PARAM.nu = 1000;
% PARAM.thresoldAM = .05;
% PARAM.thrDecay = 2^-11; % optimal around 2^-11
PARAM.learningRate_i = 2^-6; % avoid more than 2^-4
PARAM.learningRate_f = 2^-0; % avoid more than 2^0
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
%     thresholdAM = .00*ones([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat
    inputThreshold = 0.055*ones([PARAM.RFSize PARAM.nFeat]);% format: i x j x feat
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
    
    % init evol
    evol = zeros(1,PARAM.n,floor(PARAM.nIter/100));
end
disp('Learning invariances...')
tic
offset = floor(rand*length(PARAM.nIter));
% disp(int2str(offset))

for i=1:PARAM.nIter
    
    effective = false;

    % CURRENT INPUT
%     x = getMap(i+offset); % format: i x j x feat
    x = feval(getMap,i+offset); % format: i x j x feat

%     AM = x ./ trAM; % format: i x j x feat
%         AM = x; % format: i x j x feat
        input = x ./ inputThreshold; % format: i x j x feat

    [inputMax inputWinner] = max(input(:));
    inputWinnerF = floor((inputWinner-1)/PARAM.RFSize(1)/PARAM.RFSize(2))+1;
    inputWinnerJ = floor((inputWinner-(inputWinnerF-1)*PARAM.RFSize(1)*PARAM.RFSize(2)-1)/PARAM.RFSize(1))+1;
    inputWinnerI = inputWinner - (inputWinnerF-1)*PARAM.RFSize(1)*PARAM.RFSize(2) - (inputWinnerJ-1)*PARAM.RFSize(1);

    %     % Spratling normalization
    %     normFactor = sum(weight,4);
    %     for n=1:PARAM.n
    %         %                 weight(:,:,:,n) = weight(:,:,:,n) ./ sum(sum(sum(weight(:,:,:,n))));
    %         weight(:,:,:,n) = weight(:,:,:,n) ./ normFactor;
    %     end
    %     % Spratling : compute Z
    %     nodeMax = [ max(max(max(weight,[],1),[],2),[],3) ];
    %     nodeMax = nodeMax(:);
    %     inputMax = max(weight,[],4);
    %     inputMax = inputMax + ~inputMax; % avoid dividing by 0. anyway inputMax(i,j,f)=0 => for all nodes n weight(i,j,f,n)=0 => for all nodes Z(i,j,f,n)=0
    %     for n=1:PARAM.n
    %         Z = x .* ( weight(:,:,:,n).^2 ) ./ inputMax / nodeMax(n); % format: i x j x feat
    %         AT(n) = max(max(max( Z ))) / trAT(n);
    %     end

    if CORREL_PREV_INPUT_CURR_OUTPUT % correlation between previous input and current output (like Einhauser 2002)
        % compute current ouput
        for n=1:PARAM.n
%             AT(n) = max(max(max( weight(:,:,:,n).*AM ))) / trAT(n);
                    output(n) = max(max(max( weight(:,:,:,n).*input )));
        end
        % find winner in top layer
        [outputMax outputWinner] = max(output);

        % learning
        if i>1 && previousInputMax > 1%inputThreshold(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF)
            effective = true;
            % update thr
            inputThreshold(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF) = ...
                previousInputMax * inputThreshold(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF);
            % update weights:
            % store winner-winner weight
            winWeight = weight(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF,outputWinner);
            % depress everyone
            weight(:,:,:,outputWinner) = weight(:,:,:,outputWinner) + PARAM.a_neg*weight(:,:,:,outputWinner).*(1-weight(:,:,:,outputWinner));
            % potentiate synapse between winners
            weight(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF,outputWinner) = winWeight + PARAM.a_pos*winWeight*(1-winWeight);

            % reporting
            nFiringInput(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF) = ...
                nFiringInput(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF)+ 1;
            nFiringOutput(outputWinner) = nFiringOutput(outputWinner) + 1;
        end % learning

    else % correlation between previous output and current input

        % learning
        if i>1 && inputMax > 1%inputThreshold(winnerMI,winnerMJ,winnerMF)
            effective = true;
            % update thr
            %             inputThreshold(winnerMI,winnerMJ,winnerMF) = inputMax;
            inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF) = inputMax * inputThreshold(inputWinnerI,inputWinnerJ,inputWinnerF);
            % update weights:
%             if sum( [inputWinnerI inputWinnerJ inputWinnerF previousOutputWinner] ~= previousWinnerPair ) ~= 0 % different pair of winners
                % store winner-winner weight
                winWeight = weight(inputWinnerI,inputWinnerJ,inputWinnerF,previousOutputWinner);
                % depress everyone
                weight(:,:,:,previousOutputWinner) = weight(:,:,:,previousOutputWinner) + PARAM.a_neg*weight(:,:,:,previousOutputWinner).*(1-weight(:,:,:,previousOutputWinner));
                % potentiate synapse between winners
                weight(inputWinnerI,inputWinnerJ,inputWinnerF,previousOutputWinner) = winWeight + PARAM.a_pos*winWeight*(1-winWeight);

                %         % depress everyone
                %         weight(:,:,:,previousWinnerT) = (1-PARAM.learningRate) * weight(:,:,:,previousWinnerT);
                %         % potentiate synapse between winners
                %         weight(winnerMI,winnerMJ,winnerMF,previousWinnerT) =  weight(winnerMI,winnerMJ,winnerMF,previousWinnerT) + PARAM.learningRate;
                %         % normalize
                %         weight(:,:,:,winnerT) = weight(:,:,:,winnerT) / sum(sum(sum(weight(:,:,:,winnerT))));

                weight = max(weight,0); % just in case
                weight = min(weight,1); % just in case

                %                 % Spratling normalization
                %                 normFactor = sum(weight,4);
                %                 for n=1:PARAM.n
                %                     %                 weight(:,:,:,n) = weight(:,:,:,n) ./ sum(sum(sum(weight(:,:,:,n))));
                %                     weight(:,:,:,n) = weight(:,:,:,n) ./ normFactor;
                %                 end
                %
                % reporting
                nFiringInput(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF) = ... 
                    nFiringInput(previousInputWinnerI,previousInputWinnerJ,previousInputWinnerF)+ 1;
                nFiringOutput(outputWinner) = nFiringOutput(outputWinner) + 1;
                %         if sum(nFiringInput(:)>0) == length(nFiringInput(:))
                %             error('all input have fired at least once')
                %         end

                % store winner for future comparison
                previousWinnerPair = [inputWinnerI inputWinnerJ inputWinnerF previousOutputWinner];
%             end % different winner
        end % learning

        % compute current ouput
        for n=1:PARAM.n
%             AT(n) = max(max(max( weight(:,:,:,n).*AM ))) / trAT(n);
                    output(n) = max(max(max( weight(:,:,:,n).*input )));
        end
        % find winner in top layer
        [outputMax outputWinner] = max(output);

    end % switch between the 2 possible correlations

    if mod(i,100)==0
%         evol(:,:,i/100) = reshape(sum(sum(weight,1),2),[PARAM.nFeat PARAM.n]); % format: feat x node x iter
        evol(1,:,i/100) = sum(sum(sum(weight,1),2),3); % format: 1 x node x iter
    end

    % update traces
%     trAM = AM / PARAM.nu + (1-1/PARAM.nu) * trAM;
%     if effective
%         trAT = AT / PARAM.nu + (1-1/PARAM.nu) * trAT;
%     end

    % thr decay
    inputThreshold = inputThreshold*(1-PARAM.thrDecay);
    % previous := current
    previousInputMax = inputMax;
    previousInputWinnerI = inputWinnerI;
    previousInputWinnerJ = inputWinnerJ;
    previousInputWinnerF = inputWinnerF;
    previousOutputWinner = outputWinner;
    previousOutputMax = outputMax;

    if mod(i,1000)==0
        PARAM.a_pos = PARAM.a_pos * geomFactor;
        PARAM.a_neg = PARAM.a_neg * geomFactor;
        if mod(i,10000)==0
            fprintf(1,'.');
        end
    end
end
toc
% fprintf(1,'\n');

function weight = getInitialWeight(patchSize,nFeat)
% weight = 1 + .01*rand([patchSize nFeat]);
% weight = max(0,weight);
% weight = weight / sum(weight(:));
% weight = .5*( 1 + .000 * (rand([patchSize nFeat])-.5) );
weight = .5 * ones([patchSize nFeat]);

% function input = tmp(patchSize,nFeat)
% input = rand(patchSize,patchSize,nFeat);
