function ySmooth = smoothCurve(x,y, hidden, training)
%% Smooths y as function of x using a simple neural network
% Inputs
% x, y: are vectors of same length.
% hidden: (optional, default: 3) number of hidden units to be used in the neural network.
% training: (optional, default: 'trainbr') the trianing algorithm used in the neural network.
if nargin <=2
    hidden = 3;
end
if nargin <= 3
    training = 'trainbr';
end

net = feedforwardnet(hidden, training);
net.trainParam.show = NaN;
net.trainParam.showWindow = false;

net = train(net,x,y);
ySmooth = net(x);

end

