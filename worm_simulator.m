%%%%%%%
% Code to simulate tangling and untangling worms

% Plotting tool obtained from Matlab file exchange, used under the following license:
% 
% Copyright (c) 2016, Janus H. Wesenberg
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



clear all;


%Kirchhoff equation parameters
inpElastic = readmatrix('input_elastic_parameters.txt', 'range', 'A1:A9');

h = inpElastic(1); %radius of the rod
E = inpElastic(2); %elastic modulus
rho = inpElastic(3);  %density
nu = inpElastic(4); %Poisson's ratio
eta = inpElastic(5);  %2nd derivative damping parameter
zeta = inpElastic(6); %friction damping parameter
skindrag = inpElastic(7); %skin drag
D_b = inpElastic(8); %variance of random force on rod
D_tw = inpElastic(9); %variance of random torque on rod

Eb = 0.1*E;  %bending modulus


%tangle model parameters
inpHeadForce = readmatrix('input_active_head_force_parameters.txt', 'range', 'A1:A5');

chiralityNo = inpHeadForce(1);
inputAlpha = inpHeadForce(2);
poissonScale = inpHeadForce(3);
poissonD_R = inpHeadForce(4); 
poissonD_F = inpHeadForce(5); 

forceReversalRate = inputAlpha / (2*pi*chiralityNo);


%initialize dynamics

%name of output file
fileLabel = strcat(pwd, '/results/', 'worm_sim_output');

%define initial condition
inputIC = readlines('input_initial_condition.txt');
inputIC = inputIC(1);
load(strcat(pwd, '/initial_conditions/',inputIC), 'W0', 'diskR');

nw = size(W0,2)-2; worms = size(W0,1)/3; %size and number of worms
m0w = zeros(worms,nw); %initial twist density

%Other parameters which can be changed
N=40000; %number of timesteps
framenumber=200; %number of frames to capture
dt= 10^(-4); %size of timesteps



%Elastic parameters derived from inputs

mutwist=Eb/(2*(1+nu));  %shear modulus
Kbulk = 20*E/(3*(1-2*nu)); %bulk modulus
Ar=pi*h^2; %cross-sectional area
I=0.25*pi*h^4; %moment of inertia
J=0.5*pi*h^4; %moment of twist

alpha = Eb*I; %bending energy scale
beta = mutwist*J; %twisting energy scale
gamma = E*Ar; %stretching energy scale



eta_eff = eta * Ar;  %second derivative damping
eta_tw = eta * 2*I; %twist second derivative damping

zeta_eff = zeta*Ar; %collision stick factor
zeta_tw_eff = zeta*J; %twist stick factor





%rest quantities:     
[L0w, Ls0w, vLengths_fric] = initialLengths(W0);
ds = vLengths_fric'; ds = ds(:)';
massList = rho*Ar*vLengths_fric(:)';
skindrag_eff = skindrag*Ar*vLengths_fric(:)';    
mass_eff_tw = rho*2*I*L0w;
skindrag_eff_tw = skindrag*2*I*L0w;



%Loop quantities
%Need V0, velw, accelw, thetaW0, thdotW, d1w, m_twW
%     poissonAlpha, poissonN, binormals

%initial position, velocity, acceleration
V0 = W0; 
velw = zeros(3*worms, nw+2); velW = zeros(3, worms*(nw+2));
accelw = zeros(3*worms, nw+2); accelW = zeros(3, worms*(nw+2));

%initial twist, twist density, angular velocity
thetaW0 = zeros(worms, nw+1);
m_twW = m0w; 
thdotW = zeros(worms, nw+1);
d1w = initialMaterialFrame(W0);

%initialize SDE for worm head activity
[poissonAlpha, poissonN, poissonForce] = initialActiveForce(worms, inputAlpha);
binormals = [zeros(2,worms); zeros(1,worms)+1]; %set the plane for the random force SDE


%Initialize worm position array and worm twist array
W = zeros(3*worms,nw+2,framenumber+1); W(:,:,1) = W0;
TwDispW = zeros(worms,nw,framenumber+1); TwDispW(:,:,1) = m0w;


% quantities to track over time
% 'V0', 'velw', 'velW', 'accelw', 'accelW', 'thetaW0', 'm_twW',
% 'thdotW', 'd1w', 'poissonAlpha', 'poissonN', 'poissonForce'
save(strcat(fileLabel,'.mat'), 'W', 'TwDispW', 'V0', 'velw', 'velW', 'accelw', 'accelW', 'thetaW0', 'm_twW', 'thdotW', 'd1w', 'poissonAlpha', 'poissonN', 'poissonForce');



% main loop
tic; LoopTime0 = 0;
for i=1:N 


    %initialize variables, timestepping taking Y0 -> Y, etc.
    V00 = V0;
    V0 = V0 + dt*velw + 0.5*dt^2*accelw; %thetaW0 = thetaW;
    for ii = 1:worms
        
        %no penetration of plate
        V0(3*ii,V0(3*ii,:)<h) = V0(3*ii,V0(3*ii,:)<h) + 0.1*(h - V0(3*ii,V0(3*ii,:)<h));

        %cylindrical trap
        Vradial = vecnorm(V0(3*ii-2:3*ii-1,:));
        V0(3*ii-2:3*ii-1, Vradial>diskR) = (diskR./Vradial(Vradial>diskR)) .* V0(3*ii-2:3*ii-1, Vradial>diskR);

        %update twist
        thetaW0(ii,:) = thetaW0(ii,:) + dt*thdotW(ii,:);
        d1w(3*(ii-1)+1:3*ii,:) = fupdate(V00(3*(ii-1)+1:3*ii,:), V0(3*(ii-1)+1:3*ii,:), d1w(3*(ii-1)+1:3*ii,:)); %update reference frame along the rope
        m_twW(ii,:) = tAngle(V0(3*(ii-1)+1:3*ii,:), d1w(3*(ii-1)+1:3*ii,:), thetaW0(ii,:)); %new twist displacement for A

    end


    %contact handling terms
    V0comp = zeros(3, worms*(nw+2));
    for ii = 1:worms
        V0comp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = V0(3*(ii-1)+1:3*ii,:);
    end

    Rw = dists(V0comp);
    for ii = 1:worms-1 %get rid of links between worm i head and worm i+1 tail
        Rw(:,ii*(nw+2)) = 100*h; Rw(ii*(nw+2),:) = 100*h;
    end

    %change in effective radius here
    iMw = interactMat(Rw,1.25*h); pMatw = pressure(V0comp,Rw,iMw,1.25*h,Kbulk);

    %calculate elastic forces
    FelastW = zeros(3, worms*(nw+2));
    iFW = 1*iForce(V0comp,pMatw,Rw);
    for ii=1:worms
        w = V0(3*(ii-1)+1:3*ii,:);
        FelastW(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = bForce1(w,alpha) + sForce(w,L0w(ii,:),gamma) + tForce(w, m_twW(ii,:), beta);


        %construct active head force
        if exprnd(1/forceReversalRate) < dt
            poissonAlpha(ii) = -poissonAlpha(ii);
        end
        brN = poissonN(:,ii)/norm(poissonN(:,ii)); bnml = binormals(:,ii);
        d_brN = (poissonAlpha(ii) * cross(bnml, brN) - poissonD_R*brN )*dt + sqrt(dt) * sqrt(2*poissonD_R)*(eye(3) - brN*brN' - bnml*bnml')*normrnd(0,1,3,1);    
        poissonN(:,ii) = brN + d_brN;
        poissonForce(:,ii) = poissonScale * poissonN(:,ii) + sqrt(2*poissonD_F/dt)*normrnd(0,1,3,1);

        %add active head force
        FelastW(:, ii*(nw+2)) = poissonForce(:,ii) + FelastW(:, ii*(nw+2));
    end
    
    %noise term
    randForce = sqrt(2*D_b/dt)*sqrt(ds).*normrnd(0,1, 3, (nw+2)*worms);

    %total force
    FtotW = FelastW + iFW + randForce;




    %calculate hydrodynamic interaction terms implicitly in terms of
    %velocity and calculate interaction strength based on depth of penetration
    RvW = distsV(V0comp); %inverse squared distance
    iMcW = padarray(iMw, [1,1], 'post') + padarray(iMw,[1,1],'pre') - padarray(iMw,[1,1],'post').*padarray(iMw,[1,1],'pre');
    pDepthW = RvW.*iMcW; 
    pDepthW = (pDepthW .* vLengths_fric(:)) .* vLengths_fric(:);
    pDepthW = (pDepthW - diag(sum(pDepthW)));
    v0 = zeros(1,worms*(nw+2)); v1 = zeros(1,worms*(nw+2));
    for ii=1:worms
       v0((ii-1)*(nw+2)+1:ii*(nw+2)) = [1/L0w(ii,1), 1./L0w(ii,1:end-1)+1./L0w(ii,2:end), 1./L0w(ii,end)];
       v1((ii-1)*(nw+2)+1:ii*(nw+2)) = [1./L0w(ii,:), 0];
    end
    v1 = v1(1:end-1);
    Lapl0W = diag(v0) - diag(v1,1) - diag(v1,-1);
    LaplW = (eta_eff * Lapl0W - zeta_eff *pDepthW); 

    %linear operator A to get new velocity from old velocity, v_old
    A = (1+0.5*dt*skindrag/rho) * eye((nw+2)*worms) + 0.5*dt*(1./massList).*LaplW;
    v_old =  velW' + 0.5*dt*(FtotW'.*(1./massList') + accelW');

    %new velocity from solution of sparse linear system
    A = sparse(A);
    velW = lsqminnorm(A, v_old)';

    %new acceleration
    accelW = 1./massList .* (FtotW - skindrag_eff .* velW - (LaplW *velW')' );

    %timestep the positions of the ropes
    velw = zeros(3*worms,nw+2);
    accelw = zeros(3*worms,nw+2);
    for ii=1:worms
        velw(3*(ii-1)+1:3*ii,:) = velW(:,(nw+2)*(ii-1)+1:ii*(nw+2));
        accelw(3*(ii-1)+1:3*ii,:) = accelW(:,(nw+2)*(ii-1)+1:ii*(nw+2));
    end



    %calculate torques
    tauW = zeros(worms, nw+1);
    for ii = 1:worms
        tauW(ii,:) = angleForce(V0(3*(ii-1)+1:3*ii,:), m_twW(ii,:), beta);
    end

    %friction-like interaction torque
    tauIntW = zeta_tw_eff*tauForce(V0comp,velW,iMw,Rw);

    %solve implicit equation for angular velocities on each worm
    tau_totW = zeros(worms, nw+1);
    for ii = 1:worms
        w = V0(3*(ii-1)+1:3*ii,:);
        
        %noise term
        randTorque = sqrt(2*D_tw/dt) * sqrt(L0w(ii,:)) .* normrnd(0,1,1,nw+1);
        
        %total torque
        tau_totW(ii,:) = tauW(ii,:) + [0, tauIntW((ii-1)*(nw+2)+2:ii*(nw+2)-1)] + randTorque;
       
        %timestep angular velocity
        twL = eta_tw * TwLapl(w, Ls0w(ii,:));
        thdotW(ii,:) = thdotW(ii,:) + dt* 1./mass_eff_tw(ii,:) .* (tau_totW(ii,:) - skindrag_eff_tw(ii,:) .* thdotW(ii,:) - (twL*thdotW(ii,:)')' );
    end

    %save data every few timesteps
    if mod(i,N/framenumber)==0
        W(:,:,1+framenumber*i/N) = V0;
        TwDispW(:,:,1+framenumber*i/N) = m_twW;

        %save progress
        save(strcat(fileLabel,'.mat'), 'W', 'TwDispW', 'V0', 'velw', 'velW', 'accelw', 'accelW', 'thetaW0', 'm_twW', 'thdotW', 'd1w', 'poissonAlpha', 'poissonN', 'poissonForce');
        
        %estimate remaining time
        LoopTime1 = toc; LoopTime = (LoopTime1-LoopTime0) * framenumber/N * (N-i);
        fprintf('\n Estimated time remaining: %0.0f s\n', LoopTime)        
        LoopTime0 = LoopTime1;
        
    end
end
% end main loop



%make video of simulation

%%%%%create AVI object
nFrames = size(W,3);
vidObj = VideoWriter(char(strcat(fileLabel)), 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);

colormap qualitative12

%%%%create movie
for j=1:1:nFrames
    clf;
    for ii = 1:worms

        %plot worm
        [x,y,z,C]=tubeplot2(W(3*(ii-1)+1:3*ii,:,j),h,20,h/5,[zeros(1,nw+2)+mod(ii,12)]);
        Aplot=surf(x,y,z,C);
        shading interp;
        set(Aplot,'meshstyle','row');
        set(Aplot, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
        
        material shiny; 

        hold on;
        
    end
    viscircles([0,0], diskR,'Color', 'k');
    axis(diskR*[-1,1,-1,1,-1,1]);

    view(2); caxis([0,11]);
    hl=camlight('left');
    
    daspect([1,1,1]); 
    grid off; box on;
    set(gca, 'XTick',[]); set(gca, 'YTick',[]); set(gca,'ZTick',[]);
    set(gca,'Color',[0.6 0.6 0.6]);
    set(gcf, 'Color', [1,1,1]);
    set(gca,'LooseInset',get(gca,'TightInset'));

    writeVideo(vidObj, getframe(gcf));
    drawnow;
    hold off;
    fprintf('\nCompiling video %0.0f%%\n', 100*j/nFrames);

end
close(vidObj);
hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [L0w, Ls0w, vLengths_fric] = initialLengths(W0)
    nw = size(W0,2)-2; worms = size(W0,1)/3;
  
    L0w = zeros(worms, nw+1);
    Ls0w = zeros(worms, nw+2);
    for i = 1:worms
        w = W0(3*(i-1)+1:3*i,:);
        L0w(i,:) = vecnorm(w(:,2:end) - w(:,1:end-1));
        Ls0w(i,:) = Vdom(w);
    end
    
    vLengths_fric = Ls0w'; 
    vLengths_fric(1,:) = 2*vLengths_fric(1,:); vLengths_fric(end,:) = 2*vLengths_fric(end,:);
    
end


function d1w = initialMaterialFrame(W0)
    nw = size(W0,2)-2; worms = size(W0,1)/3;
    d3w = W0(:,2) - W0(:,1);
    d1w = zeros(3*worms,nw+1);
    for i = 1:worms
        w = d3w(3*(i-1)+1:3*i);
        d1 = cross(cross(w, [0,0,1]'), w); 
        if norm(d1)==0
            d1 = cross(cross(w, [0,0,1.01]'), w);
        end
        d1=d1./norm(d1);
        d1w(3*(i-1)+1:3*i,:) = pTransport(W0(3*(i-1)+1:3*i,:), d1);
    end
end


function [poissonAlpha, poissonN, poissonForce] = initialActiveForce(worms, alpha)
    poissonAlpha = alpha + zeros(1,worms);
    th = 2*pi*rand(1,worms);
    poissonN = [cos(th); sin(th); zeros(1,worms)];
    poissonForce = zeros(3, worms);
end




function R = dists(X) %distance between i'th and j'th links
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    L2 = diag(Xm'*Xm);
    R = -2*(Xm'*Xm) + (L2 + L2');
    R = R.^0.5; 
end

function R = distsV(X) %inverse squared distance between i'th and j'th points
    L2 = diag(X'*X);
    R = -2*(X'*X) + (L2 + L2');
    R = 1./(R + eye(size(R))) - eye(size(R));
end




function iM = interactMat(R, h) %Set R(i,j)=0 for |i-j|<=10
    n = size(R)-1;
    R([1:n+2:end, n+1+1:n+2:end-1]) = 2*h; 
    R([2:n+2:end-(n+1)]) = 2*h; 
    iM = R<2*h; %iM(i,j) = 1 iff link i and j are interacting
end


function pMat = pressure(X,R,iM,h,Kbulk)
    R(iM==0) = 0;
    pMat = Kbulk*( 1*(1 - R./(2*h)) + 100*(1-R./(2*h)).^3 + 100*(1-R./(2*h)).^5 + 100*(1-R./(2*h)).^7); %set cubic term to around 100
    d3 = X(:,2:end) - X(:,1:end-1);
    L = abs(d3'*d3); L(iM==0)=0;
    pMat = pMat.*L;
end





function Fi = iForce(X, pMat, R)
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    R(R==0) = 1;
    pMat = pMat./R;
    pMat = sum(pMat).*Xm - (Xm*pMat);
    pMat=[pMat(:,1), 0.5*(pMat(:,2:end)+pMat(:,1:end-1)), pMat(:,end)];
    Fi=pMat;
end



function L = TwLapl(X, Ls)
    ph = phi(X); ph = cos(ph(2:end-1));
    L = diag([1./Ls(2), 1./Ls(2:end-2)+1./Ls(3:end-1),1./Ls(end-1)]) - diag([ph(1:end)./Ls(2:end-1)],1) - diag(ph(1:end)./Ls(2:end-1),-1);
end


function tau = tauForce(X, F, iM, R)
    d3 = X(:,2:end) - X(:,1:end-1); %d3 = d3./vecnorm(d3);
    iM = iM./(eye(size(R)) + R.^4);
    iM = (iM .* vecnorm(d3)') .* vecnorm(d3); d3 = d3./vecnorm(d3);
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    Fm = 0.5*(F(:,1:end-1)+F(:,2:end));
    tau = cross(Xm, Fm)*iM - cross(Xm, Fm*iM);
    tau = tau - cross(Xm*iM, Fm) + full(sum(iM)).*(cross(Xm, Fm));
    tau = dot(tau, d3);
end


function B = pTransport(X,v0) %X is the set of positions, 3 by n+2, %v0 is initial vector
    n=size(X);
    n=n(2)-2;
    Y=zeros(3,n+1);
    Y(:,1)=v0; %initial frame
    for i=2:n+1
        ti=X(:,i)-X(:,i-1);
        tf=X(:,i+1)-X(:,i);
        R = rot(ti,tf);
        Y(:,i)=R*Y(:,i-1);
    end
    B=Y;
end



function R = rot(ti,tf)
        theta=atan2(norm(cross(ti,tf)),dot(ti,tf)); %angle to rotate ti into tf
        if norm(cross(ti,tf))==0
            m=[0,0,0];
        else
            m=cross(ti,tf)./norm(cross(ti,tf)); %the sign of m comes from here
        end
        m=theta*m;
        A=[0,-m(3),m(2);m(3),0,-m(1);-m(2),m(1),0]; %infinitesimal rotation
        R=expm(A);
end



function kb = darboux(X)
    d3=(X(:,2:end)-X(:,1:end-1));
    t=d3./sqrt(dot(d3,d3));
    kb=2*cross(t(:,1:end-1),t(:,2:end))./(1+dot(t(:,1:end-1),t(:,2:end)));
end



function Fs = sForce(X0,L0,gamma)
    d3 = X0(:,2:end)-X0(:,1:end-1);
    Ln = vecnorm(d3); d3 = d3./Ln;
    stretch = Ln./L0 - 1;
    F = gamma.*(stretch + 0*stretch.^3).*d3;
    Fs = [F, [0,0,0]'] + [[0,0,0]', -F];
end



function ph = phi(X) %turning angle from one segment to the next 
    n=size(X,2)-2;
    d3=X(:,2:n+2)-X(:,1:n+1);
    Y=atan2(sqrt(sum(cross(d3(:,1:n),d3(:,2:n+1)).^2,1)),dot(d3(:,1:n),d3(:,2:n+1)));
    ph=[0,Y,0];
end



function Ls = Vdom(X) %vector containing size of voronoi domain at each X
    n=size(X,2)-2;
    Y=(X-[X(:,1),X(:,1:n+1)]);
    Z=([X(:,2:n+2),X(:,n+2)]-X);
    Ls=0.5*(sqrt(sum(Y.^2))+sqrt(sum(Z.^2)));
end



        
function dp = dphi(X) %dp(i,j,:) = dphi(i)/dx(j)
    n=size(X);
    n=n(2)-2;
    Yp=zeros(n+2,n+2,3);
    d3=X(:,2:n+2)-X(:,1:n+1);
    cosph = (dot(d3(:,1:n), d3(:,2:n+1))./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))));
    dcos = -1./sqrt(1-cosph.^2);
    ld = -d3(:,2:n+1)./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))) + (d3(:,1:n)./(sum(d3(:,1:n).^2))).*cosph;
    ud = d3(:,1:n)./(sqrt(sum(d3(:,1:n).^2)).*sqrt(sum(d3(:,2:n+1).^2))) - (d3(:,2:n+1)./(sum(d3(:,2:n+1).^2))).*cosph;
    d=[[0,0,0]',-dcos.*(ld+ud),[0,0,0]'];
    ld=[dcos.*ld,[0,0,0]']; ud=[[0,0,0]',dcos.*ud];
    Yp(:,:,1)=(diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1));
    Yp(:,:,2)=(diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1));
    Yp(:,:,3)=(diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1));
    Yp=real(Yp);
    Yp(isnan(Yp))=0; Yp(isinf(Yp))=0;
    dp=Yp;
end


function derivs = dL(X) %derivs(i,j,:) = dL(i)/dx(j) where L(i) is twice the voronoi domain of i'th vertex
    n=size(X);
    n=n(2)-2;
    Yp=zeros(n+2,n+2,3);
    ld=(X(:,1:n+1)-X(:,2:n+2))./(sqrt(sum((X(:,1:n+1)-X(:,2:n+2)).^2)));
    ud=-ld;
    d=[[0,0,0]',ud]-[ud,[0,0,0]'];
    Yp(:,:,1)=diag(d(1,:))+diag(ud(1,:),1)+diag(ld(1,:),-1);
    Yp(:,:,2)=diag(d(2,:))+diag(ud(2,:),1)+diag(ld(2,:),-1); 
    Yp(:,:,3)=diag(d(3,:))+diag(ud(3,:),1)+diag(ld(3,:),-1);
    derivs=Yp;
end



function Fb=bForce1(X,alpha)
    n=size(X);
    n=n(2)-2;
    ph=phi(X); dph=dphi(X);
    L=2*Vdom(X); dl=dL(X);
    Yb=zeros(3,n+2);
    Yb(1,:)= -((2*alpha*ph)./L)*dph(:,:,1) + ((alpha*ph.^2)./L.^2)*dl(:,:,1);
    Yb(2,:)= -((2*alpha*ph)./L)*dph(:,:,2) + ((alpha*ph.^2)./L.^2)*dl(:,:,2);
    Yb(3,:)= -((2*alpha*ph)./L)*dph(:,:,3) + ((alpha*ph.^2)./L.^2)*dl(:,:,3);
    Fb=Yb;
end

    
function m = tAngle(X, d1, theta) %d1(i) is reference frame vector, theta(i) is angle on link i 
    d3=X(:,2:end)-X(:,1:end-1);
    N=cross(d3(:,1:end-1),d3(:,2:end));
    N=N./sqrt(dot(N,N));
    N(isnan(N))=0;
    d3=d3./sqrt(dot(d3,d3));
    ph=phi(X);
    ph=ph(2:end-1);
    pTrans=cos(ph).*d1(:,1:end-1)+((1-cos(ph)).*dot(N,d1(:,1:end-1))).*N + (sin(ph).*cross(N,d1(:,1:end-1)));  
    pTrans(isnan(pTrans))=d1(isnan(pTrans));
    frameTwist=atan2(dot(d1(:,2:end),cross(d3(:,2:end),pTrans)),dot(d1(:,2:end),pTrans));
    m=theta(2:end)-theta(1:end-1)+frameTwist;
end



function Ft = tForce(X,m,beta) %force on x_i due to twisting, m is given by tAngle
    Ls=Vdom(X);
    kb=darboux(X);
    Es=sqrt(dot(X(:,2:end)-X(:,1:end-1),X(:,2:end)-X(:,1:end-1))); %Es(i) is distance between x_i and x_i+1
    %forces on x_i due to changing the frame due to changing m at i-1, i+1,
    %i respectively
    Fd=0.5*[[0,0,0]', [0,0,0]', (m(1:end)./Ls(2:end-1)).*(kb(:,1:end)./Es(2:end))]; 
    Fu=0.5*[(m(1:end)./Ls(2:end-1)).*(-kb(:,1:end)./Es(1:end-1)), [0,0,0]',[0,0,0]'];
    Fe=0.5*[[0,0,0]', (m(1:end)./Ls(2:end-1)).*(kb(:,1:end)./Es(1:end-1) - kb(:,1:end)./Es(2:end))  ,[0,0,0]'];
    twist1=-beta*(Fd+Fu+Fe);
    %force on x_i due to changing the lengths of links in the twist energy
    %expression
    L=2*Ls; dl=dL(X); 
    twist2=[((beta*[0,m,0].^2)./L.^2)*dl(:,:,1);((beta*[0,m,0].^2)./L.^2)*dl(:,:,2);((beta*[0,m,0].^2)./L.^2)*dl(:,:,3)];
    Yt=twist1+twist2;
    Ft=Yt;
end


function Fa=angleForce(X,m,beta)
    Ls=Vdom(X);
    Fa= -beta*([0,m./Ls(2:end-1)]-[m./Ls(2:end-1),0]);
end

function d1New = fupdate(X0,X1,d1)
    d30=X0(:,2:end)-X0(:,1:end-1); d31=X1(:,2:end)-X1(:,1:end-1);
    ph=atan2(sqrt(sum(cross(d30,d31).^2,1)),dot(d30,d31));
    N=cross(d30,d31);
    N=N./sqrt(dot(N,N));
    Y=cos(ph).*d1+((1-cos(ph)).*dot(N,d1)).*N + (sin(ph).*cross(N,d1));  
    Y(isnan(Y))=d1(isnan(Y));
    d1New=Y;
end




function [x,y,z,C]=tubeplot2(curve,r,n,ct,S)
% Usage: same as above but this gives colours the rod according to a 
% scalar field S along the curve.

  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
     end
  end
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end
  if nargin<4 || isempty(ct)
    ct=0.5*r;
  end

  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end

  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [~,idx]=min(abs(dv(:,1))); nvec(idx)=1;

  xyz=zeros(3,n+1,npoints+2);
  Col=zeros(3,n+1,npoints+2); 

  %precalculate cos and sing factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
    Col(:,:,k+1)=S(k);
  end
  %finally, cap the ends:
  xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  Ct=squeeze(Col(1,:,:));
  C=Ct; C(:,1) = C(:,2); C(:,end) = C(:,end-1);
  %... and plot:
  if nargout<3, surf(x,y,z,C); end
end





