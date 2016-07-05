function [ Gz ] = RecursiveLeastSquarewithVaryingExponentialForgetting ( u, y, na, nb, d, Ts, Plot )
% this function is used to estimate the parameter of the system without
% noise by using RLS method with time varying parameter by using factor lambda
% estimated G(z)=z^d(b_0 z^nb+b_1 z^(nb-1)+...+b_nb)/(z^na+a_1 z^(na-1)+...+a_na)
% Y=phi'*theta_hat
% where :
% phi'=[y(t-1) ... y(t-na) u(t-1-d) ... u(t-nb-d)]
% theta_hat =[a_1 ... a_na b_0 ... b_nb]
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs:
% uc    : input required signal
% y      : output signal to the system
% na   : order of thr Den. of the transfer function of the system
% ab   : order of thr Num. of the transfer function of the system
% d     : order of delay of the system
% Ts     : sampling time
% Plot   :  if Plot(1) = 0 -> No plot
%             if Plot(1) = 1 ->  plot    Plot (2:4) = number of figures
%% outputs:
% Gz              : discreate transfer function
%% function body
N=max(na,nb+d+1);
P{1}=eye(na+nb+1)*10^(4);
theta_hat(:,1)=zeros(na+nb+1,1);
lambda(1)=0.9;
xxx=1;
for i=2:length(y)
    for j=1:na
        if i-j <=0
            phiT(i,j)=0;
        else
            phiT(i,j)=[-y(i-j)];
        end
    end
    for j=0:nb
        if i-j-d <= 0
            phiT(i,j+1+na)=0;
        else
            phiT(i,j+1+na)=[u(i-j-d)];
        end
    end
    K1{i}=P{i-1}*phiT(i,:)'*inv(lambda(i-1)+phiT(i,:)*P{i-1}*phiT(i,:)');
    epslon(i)=y(i)-phiT(i,:)*theta_hat(:,i-1);
    theta_hat(:,i)=theta_hat(:,i-1)+K1{i}*epslon(i);
    P{i}=(eye(na+nb+1)*lambda(i-1)-K1{i}*phiT(i,:))*P{i-1};
    if i > N
            lambda(i)=1-(1-phiT(i,:)*K1{i})*epslon(i)^2/(std(epslon)^2*mean(epslon));
    else
            lambda(i)=lambda(i-1);
    end
end

Theta_hat=theta_hat(:,end);
Gz=tf([Theta_hat(na+1:end)'],[1,Theta_hat(1:na)'],Ts);

for l=1:length(phiT(:,1))
    y1(l)=phiT(l,:)*theta_hat(:,l);
end
%% Plotting
if Plot(1)==1
    figure(Plot(2));
    set(gcf,'color','w')
    set(0,'defaultfigureposition',[40 55 1300 600])
    hold all;
    for k=1:na
        plot((0:length(y1)-1)*Ts,theta_hat(k,:),'linewidth',2);
    end
    for k=1+na:na+nb+1
        plot((0:length(y1)-1)*Ts,theta_hat(k,:),'--','linewidth',2);
    end
    grid on;
    for m=1:na
        Ylabel{m}=['a_' num2str(m) ', '];
        Leg{m}=['a_' num2str(m) ];
    end
    for n=1:nb+1
        if n<nb+1
            Ylabel{n+na}=['b_' num2str(n-1) ', '];
            Leg{n+na}=['b_' num2str(n-1) ];
        else
            Ylabel{n+na}=['b_' num2str(n-1)];
            Leg{n+na}=['b_' num2str(n-1) ];
        end
    end
    xlabel('t(s)','fontsize',18);
    Ylabel=cell2mat(Ylabel);
    ylabel(Ylabel,'fontsize',18);
    legend(Leg)
    title('Parameter estimation with time','fontsize',18)
    
    figure(Plot(3));
    subplot(3,1,1:2)
    set(gcf,'color','w')
    plot((0:length(y1)-1)*Ts,y,(0:length(y1)-1)*Ts,y1,'-o','linewidth',2);
    grid on;
    ylabel('y, y_e_s_t','fontsize',18);
    legend('y','y_e_s_t')
    subplot(3,1,3)
    plot((0:length(y1)-1)*Ts,abs(y-y1))
    max([abs(y-y1)]);
    grid on;
    xlabel('t(s)','fontsize',18);
    ylabel('y-y_e_s_t','fontsize',18);
    legend('error')
    
    figure(Plot(4));
    set(gcf,'color','w')
    plot((0:length(lambda)-1)*Ts,lambda,'linewidth',2);
    grid on;
    ylabel('\lambda (t)','fontsize',18);
    xlabel('t(s)','fontsize',18);
end
end