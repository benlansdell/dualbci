close all
clear
clc



k_len=10; % length of the input filters
k_mat(1,:)=.1*exp(2*(1:k_len)/k_len).*cos(2*(1:k_len)/k_len*2*pi); % arbitrarily chosen input filters
k_mat(2,:)=.05+.1*cos(2*pi*(1:k_len)/k_len);

h=[-4 -3 -2 -0.5 0.5 1 1 0.5 -0.5 -0.5]; % arbitrarily chosen feedback filter
h_len=length(h); % h is a row vector, while k is a column vector...why? no reason.  

figure
subplot(2,1,1)
plot(k_mat')
subplot(2,1,2)
plot(h)

b=-6; % basal lever pressing rate = logistic(b)


n_steps=20000; % number of time steps, maybe dt=50 ms


x_mat(1,:)=ceil(rand(1,n_steps)*5).*(rand(1,n_steps)<.05); % input matrix, could be binned spike rates, or LFP power
x_mat(2,:)=ceil(rand(1,n_steps)*2).*(rand(1,n_steps)<.05); % this could be the second cell or LFP channel
                                                           % sparse random integers in this case
                                                           % each row is a channel, and each column is a time step

[n_input x_len]=size(x_mat); % get the number of input channels and the number of time steps.  x_len==n_steps


rho=zeros(1,x_len); % preallocate the output lever press probability 
y=zeros(1,x_len); % preallocate the output lever presses

t_vec=(1:x_len); %*dt;

for i=1:x_len % step through time

    k_ind = 1:min(k_len,i); % index the k filters.  % At the beginning, the filters aren't entirely full, so min makes sure the index doesn't go negative
    x_ind_k=max(1,i-k_len+1):i; % index the input matrix.  % Some kludgy indexing to make sure the beginning and end work out

    h_ind = 1:min(h_len,(i-1)); % indexing the h filter, similar to indexing the k filter
    y_ind_h=max(1,i-h_len+1-1):(i-1); % indexing the history of the lever presses, note that it is (i-1) so that it doesn't have access to the current time step, which is defined to be in the future (i.e. the next time step).
    
    
    sum_in_nonlin=0; % sum the stuff that goes into the nonlinearity
    sum_in_nonlin=sum_in_nonlin+b; % first, add the constant lever press probability term
    for i_input=1:n_input

        if i>1 % convolve the k filters and convolve the lever press filters, add them together
            sum_in_nonlin=sum_in_nonlin+dot(fliplr(k_mat(i_input,k_ind)),x_mat(i_input,x_ind_k)) + dot(fliplr(h(h_ind)),y(y_ind_h));
        else % at i<1, there is no lever history
            sum_in_nonlin=sum_in_nonlin+dot(fliplr(k_mat(i_input,k_ind)),x_mat(i_input,x_ind_k));
        end
        
    end
    
    rho(i)=1/(1+exp(-sum_in_nonlin)); % the output after the nonlinearity
    
    y(i) = binornd(1,rho(i)); % sample the distribution for the current time step to find out if the lever was pressed
end

n_lever_presses=sum(y) % print to screen the total number of lever presses

figure 
for i_input=1:n_input
    plot(t_vec,7*(i_input-1)+x_mat(i_input,:),'b') % plot the input signals in blue
    hold on
end
plot(t_vec,rho,'r') % plot rho in red
hold on
plot(t_vec,y,'g') % plot the lever presses in green
if max(y)>30 % some times it blows up.  Then it needs to be re run.  This make the plot reasonable
    ylim([-3 12]) 
end

x=x_mat'; % flip things about before saving
y=y';
k=k_mat';
h=h';
rho=rho';
save('generated_data_1.mat','x','y','k','b','rho','h')

