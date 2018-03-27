clear all;
mu = 0.8;
N=80;
M=floor(mu*N/2);
%M=30;
r = linspace(-pi/2,pi/2,M);
theta = (0:N-1)/N * 2 * pi;
m = abs(randn(M,N));
%m(11:20,40:50)=5;
m(1:10,30:40)=5;
m_next = m;
I = 1.0;
IE = @(x) cos(x);
J1=50;
J0=-50;
W = @(t1,r1,t2,r2) J1/2*(cos(t1-t2-mu*(r1-r2)) + cos(t1-t2+mu*(r1-r2)))+J0;
[R,THETA] = ndgrid(r,theta);
theta_A = mod(THETA-mu*R,2*pi);
theta_B = mod(THETA+mu*R,2*pi);
dt=0.2;
epsilon = 0.5;
h1 = figure(1);
set(gcf, 'Position', [200, 200, 400, 300])
h2 = figure(2);
set(gcf, 'Position', [800, 200, 200, 400])
filename1 = 'cylinder_input.gif';
filename2 = 'cylinder_input2.gif';
show_original_frame = false;
num = numel(m);
shuf = randperm(num);
for t=0:dt:15
    t
    min(min(m))
    max(max(m))
    mean(mean(m))
    
    %Psi_t = min(pi+2*pi*t/10,2*pi+pi);
    Psi_t = pi;
    It = (I+epsilon)*IE(THETA-mu*R-Psi_t);
    
    %if ~show_original_frame
    figure(1);
    subplot(2,1,1)
    imshow(m,[],'Colormap',jet,'InitialMagnification',300);
    axis image manual xy
    axis off
    subplot(2,1,2)
    imshow(It,[],'Colormap',jet,'InitialMagnification',300);
    axis image manual xy
    axis off
    %else
    figure(2);
    subplot(2,1,1)
    colormap jet
    scatter(theta_A(shuf),theta_B(shuf),10,m(shuf),'filled','diamond');
    axis([0, 2*pi, 0, 2*pi])
    axis image manual
    axis off
    subplot(2,1,2)
    scatter(theta_A(shuf),theta_B(shuf),10,It(shuf),'filled','diamond');
    axis([0, 2*pi, 0, 2*pi])
    axis image manual
    axis off
    %end
    if true
        drawnow
        frame1 = getframe(h1);
        frame2 = getframe(h2);
        im1 = frame2im(frame1);
        im2 = frame2im(frame2);
        [imind1,cm1] = rgb2ind(im1,256);
        [imind2,cm2] = rgb2ind(im2,256);
        if t == 0
            imwrite(imind1,cm1,filename1,'gif', 'Loopcount',inf,'DelayTime',0.1);
            imwrite(imind2,cm2,filename2,'gif', 'Loopcount',inf,'DelayTime',0.1);
        else
            imwrite(imind1,cm1,filename1,'gif','WriteMode','append','DelayTime',0.1);
            imwrite(imind2,cm2,filename2,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
    
    
    for i=1:M
        for j=1:N
            t1 = theta(j);
            r1 = r(i);
            w = W(t1,r1,THETA,R);
            s = sum(w(:).*m(:))/N/M ;
            deriv = -m(i,j)+max(0,s+It(i,j));
            m_next(i,j) = m(i,j) + dt * deriv;
        end
    end
    m=m_next;
end