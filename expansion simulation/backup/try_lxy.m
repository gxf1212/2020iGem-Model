% ��ɫ����
Map = [1 1 1; 0 0 0];
colormap(Map);
% ���������С
S = 121;
L = zeros(S);
K = eye(S);
% ���м�һ��������Ϊ 1 ��ΪԪ������
M = (S+1)/2;
L(M, M) = 1;
Temp = L;
imagesc(L);
% �������
Layer = (S-1)/2 + 1;
%�����ʺ�������
alpha=0.1;beta=0.05;gamma=0.05;
f=getframe(gcf);
imind=frame2im(f);
[imind,cm]=rgb2ind(imind,256);
imwrite(imind,cm,'1.gif','gif','Loopcount',inf,'DelayTime',0.1);
for t=1:Layer
    for x=M-t+1:M+t-1
          for y=M-t+1:M+t-1
              if  Temp(x,y)==1 
                  i=1; 
                  A=rand(10);
                  for m=-1:1
                      for n=-1:1
                          g=A(i);
                           %g����growth�������������ģ�ֳⷱ����
                  if g <= alpha && g>=0 && Temp(x+m,y+n)==0 && K(x+m,y+n)==1 
                      %����һ��ɹ���ɢ
                               Temp(x+m,y+n)=1;
                  end
                  i=i+1;
                  end
                  end
                              K(x,y)=K(x,y)-0.015;
                             if K(x,y)<=0
                                 %�˴���Դ�ѱ�������
                                 K(x,y)=0;
                               if K(x,y)==0;
                                 nd=rand(1);
                             %nd����die of nutrition�������������ģ������Դȱ���������ĸ���
                             if nd <= gamma && nd>=0
                                 %�˴�����������
                                 Temp(x,y)=0;
                             end
                               end
                             d=rand(1);
                             %d����death�������������ģ����������
                             if d <= beta && d>=0
                                 %�˴�����������
                                 Temp(x,y)=0;
                             end
                             end
              end
          end
    end
    L = Temp;
  imagesc(L);
    % �ٶȿ���
   pause(0.5);
   %����gifͼ
   f=getframe(gcf);
   imind=frame2im(f);
   [imind,cm]=rgb2ind(imind,256);
   imwrite(imind,cm,'1.gif','gif','WriteMode','append','DelayTime',0.1);
end
