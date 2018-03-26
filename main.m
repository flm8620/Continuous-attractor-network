clear all;
mu = 0.8;
N=80;
%M=floor(mu*N/2);
M=30;
r = linspace(-pi/2,pi/2,M);
theta = (0:N-1)/N * 2 * pi;
m = abs(randn(M,N));
%m(10:20,1:end)=1;
m_next = m;
I = 1.0;
IE = @(x) cos(x);
J1=50;
J0=-50;
W = @(t1,r1,t2,r2) J1/2*(cos(t1-t2-mu*(r1-r2)) + cos(t1-t2+mu*(r1-r2)))+J0;
[R,THETA] = ndgrid(r,theta);
theta_A = mod(THETA-mu*R,2*pi);
theta_B = mod(THETA+mu*R,2*pi);
dt=0.1;
epsilon = 1;
h = figure(1);
set(gcf, 'Position', [200, 200, 400, 200])
figure(2);
set(gcf, 'Position', [800, 200, 200, 200])
axis tight manual
filename = 'moving_single_ring_stop.gif';
show_original_frame = false;
num = numel(m);
shuf = randperm(num);
for t=0:dt:100
    t
    min(min(m))
    max(max(m))
    mean(mean(m))
    %if ~show_original_frame
        figure(1);
        imshow(m,[],'Colormap',jet,'InitialMagnification',300);
        axis image manual
        axis off
    %else
        figure(2);
        colormap jet
        scatter(theta_A(shuf),theta_B(shuf),10,m(shuf),'filled','diamond');
        axis([0, 2*pi, 0, 2*pi])
        axis image manual
        axis off
    %end
    drawnow 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if t == 0 
        %imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
    else 
        %imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
    end 
    
    Psi_t = min(pi+2*pi*t/50,2*pi+pi); %pi*3/2; %min(2*pi+pi*3/2,2*pi*t/10);
    
    It = (I+epsilon)*IE(THETA-mu*R-Psi_t);
    for i=1:M
        for j=1:N
            t1 = theta(j);
            r1 = r(i);
            w = W(t1,r1,THETA,R);
            s = sum(w(:).*m(:))/N/M ;
            deriv = -m(i,j)+max(0,s+I);
            m_next(i,j) = m(i,j) + dt * deriv;
        end
    end
    m=m_next;
end