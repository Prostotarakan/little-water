clear all
V=VideoWriter('cap15rain.mp4','MPEG-4');  %сохраняет в ту же папку, где код
%V=VideoWriter('cap13rain.avi');
open(V);
S=200 ; %определяет размер изображения
x=-S:1:S;
y=-S:1:S;
[X Y]=meshgrid(x,y);
dt=0.5;%0.5;  %дискретность времени
maxt=100;  %время видео
N=0.85; %коэф.затухания
%N=N*dt;
Kk=10;  %количество капель
for k=1:Kk
    kx(k)=randi([-S+10,S-10],1); %положение капли
    ky(k)=randi([-S+10,S-10],1);
    kz(k)=-rand(1)*2-1; %начальная амплитуда
    v(k)=rand(1)*2+3; %скорость капли
    T(k)=rand(1)*maxt*3/4; %время падения капли
    t=1; %действительное время
    om(k)=v(k)/2/pi;% частота колебания
    A(:,:,k)=zeros(length(x),length(y));
    r(:,:,k)=ones(length(x),length(y))*10*S;
    fi(:,:,k)=zeros(length(x),length(y));
end;
T(k)=1;
for t=0:dt:maxt
    for k=1:Kk 
        if t>T(k)&& mod(t,1)==0
            kz(k)=kz(k)*N;
            A(k)=A(k)*N;
        end
        for i=1:length(x) 
            for j=1:length(y)
                R(i,j,k)=sqrt((x(i)-kx(k))^2+(y(j)-ky(k))^2); %Расчет радиуса вообще
                if (R(i,j,k)<(t-T(k))*v(k))&&(r(i,j,k)==10*S) 
                    if R(i,j,k)<10
                        r(i,j,k)=10;
                    else
                        r(i,j,k)=R(i,j,k); %радиус действия
                    end
                    %ваще гуд%%A(i,j,k)=kz(k)/r(i,j,k); %амплитуда в нем
                    %%%A(i,j,k)=kz(k)*exp(-r(i,j,k)/10);%!%%kz(k)/((r(i,j,k)^(N/2))+1); %САВСЭМ ХАРАШО%%амплитуда в нем
                    %A(i,j,k)=kz(k)*0.45*((r(i,j,k)+exp(1))^(1/(r(i,j,k)+exp(1))));
                    A(i,j,k)=kz(k)*((r(i,j,k)+exp(1))^(1/2/(r(i,j,k)+exp(1))));
                    fi(i,j,k)=t; %запоминаем его сдвиг
                end
            end;
        end;
    end;
    Z=zeros(length(x),length(y));
    for h=1:Kk
        Z=Z+A(:,:,h).*sin(om(h).*(t-fi(:,:,h)));
    end
    Buf=figure(1); 
    surfl(X,Y,Z)
    shading interp;
    zlim([-3,6])
    grid off
    axis ('manual');
    colormap bone
    frame=getframe(Buf);
    writeVideo(V,frame);
    
end
close(V)