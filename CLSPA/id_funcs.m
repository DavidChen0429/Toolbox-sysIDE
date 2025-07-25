function [Ga,sys] = id_funcs(u,y,dt,resamp,n)

[Ga,ws, Coh, mSuu,mSyy,mFuu, mFyy] = spa_avf(u,y,dt,10,[],[],'hamming');
Ga = frd(Ga,ws);
[HHm,HHp,wss]=bode(Ga); wss=wss/2/pi;
%%
% Defining a number of constants
p = 200;     % past window size
f = 200;     % future window size
%resamp=20;
ud3=u(resamp:resamp:end,:);
yd3=y(resamp:resamp:end,:);
%yd3=lowpass(yd3p,5,1/dt);
% PBSID-opt
[S,x] = dordvarx(ud3,yd3,f,p,'tikh','gcv');
figure, semilogy(S,'x','MarkerSize', 10);
title('Singular values')
%n=6
x = dmodx(x,n);
[Ai,Bi,Ci,Di] = dx2abcdk(x,ud3,yd3,f,p,'stable1');
%%
sys=ss(Ai,Bi,Ci,Di,dt*resamp);
x0 = dinit(Ai,Bi,Ci,Di,ud3,yd3);
[HHm,HHp,wss]=bode(sys); wss=wss/2/pi;
%ysim=lsim(sys,ud3,[]); vaf((ysim),(detrend(yd3)))

% % %% moesp
% [Sn,Rnew]=dordpo(ud3,yd3,p);
% [Ae,Ce,Ke]=dmodpo(Rnew,n,'stable');
% [Be,De]=dac2bd(Ae,Ce,ud3,yd3);
% sysm=ss(Ae,Be,Ce,De,dt*resamp);
% x0 = dinit(Ae,Be,Ce,De,ud3,yd3);
% ysim2=lsim(sysm,ud3,[],x0); vaf(ysim2,yd3)

%figure;plot(ysim); hold on;plot(detrend(yd3));%plot(detrend(ysim2))
%figure;pwelch((ysim),[],[],[],1/.0005/resamp);hold on;pwelch(detrend(yd3),[],[],[],1/dt/resamp)

end