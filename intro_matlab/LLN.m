%%%%LLN%%%%%%%%

number_of_trials=10:10:10000;
%linear space that starts from 10 to 10000 but with width 10

x_bar=[];

for j=1:length(number_of_trials)
    rng(1);% seed number    
    unif_j=rand(number_of_trials(j),1);
    %rand: randomly generated number from a unifrom distribution
    %randn: randomly generated number from a standard normal distribution
    toss_j=(unif_j<=0.5);
    x_bar(j)=mean(toss_j);
end

true=0.5*ones(1,length(number_of_trials));
plot(number_of_trials,x_bar,number_of_trials, true)
title('Law of Large Numbers')
legend('sample mean','population mean')