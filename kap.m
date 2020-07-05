clear all
V=VideoWriter('cap15rain.mp4','MPEG-4');  %��������� � �� �� �����, ��� ���
%V=VideoWriter('cap13rain.avi');
open(V);
S=200 ; %���������� ������ �����������
x=-S:1:S;
y=-S:1:S;
[X Y]=meshgrid(x,y);
dt=0.5;%0.5;  %������������ �������
maxt=100;  %����� �����
N=0.85; %����.���������
%N=N*dt;
Kk=10;  %���������� ������
for k=1:Kk
    kx(k)=randi([-S+10,S-10],1); %��������� �����
    ky(k)=randi([-S+10,S-10],1);
    kz(k)=-rand(1)*2-1; %��������� ���������
    v(k)=rand(1)*2+3; %�������� �����
    T(k)=rand(1)*maxt*3/4; %����� ������� �����
    t=1; %�������������� �����
    om(k)=v(k)/2/pi;% ������� ���������
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
                R(i,j,k)=sqrt((x(i)-kx(k))^2+(y(j)-ky(k))^2); %������ ������� ������
                if (R(i,j,k)<(t-T(k))*v(k))&&(r(i,j,k)==10*S) 
                    if R(i,j,k)<10
                        r(i,j,k)=10;
                    else
                        r(i,j,k)=R(i,j,k); %������ ��������
                    end
                    %���� ���%%A(i,j,k)=kz(k)/r(i,j,k); %��������� � ���
                    %%%A(i,j,k)=kz(k)*exp(-r(i,j,k)/10);%!%%kz(k)/((r(i,j,k)^(N/2))+1); %������ ������%%��������� � ���
                    %A(i,j,k)=kz(k)*0.45*((r(i,j,k)+exp(1))^(1/(r(i,j,k)+exp(1))));
                    A(i,j,k)=kz(k)*((r(i,j,k)+exp(1))^(1/2/(r(i,j,k)+exp(1))));
                    fi(i,j,k)=t; %���������� ��� �����
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