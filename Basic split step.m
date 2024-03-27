load("examenf3.mat")
close all; clc;
N=1024; L=40; dx=L/N; dk=2*pi/L;
x=t*dx; 
k=[-N/2:1:N/2-1]*dk;
kshift=fftshift(k);kshift2=kshift.^2;
dz=dx.^2/4; zfinal=2;
pasos=ceil(zfinal/dz);

uo=u;



utotal=zeros(N,pasos+1);utotal(:,1)=uo;
un=uo;

for cuenta=1:1:pasos

  F_NL=fft(exp(1i*dz*abs(un).^2).*un);
  F_D= exp(-1i*kshift2*dz/2).*F_NL;
  un=ifft(F_D);
  utotal(:,cuenta)=un;
    
end
distancia = [0:pasos]*dz;

figure(2);
imagesc(distancia,x,abs(utotal)); colorbar;
xlabel('z, distancia'); ylabel('t, tiempo');

figure(1);
plot(x,abs(uo),x,abs(un),'r');
legend('Pulso inicial','Pulso final');
ylabel('Amplitud'); xlabel('tiempo');